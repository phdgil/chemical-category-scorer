from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from rdkit.Chem import QED

from build_g006_final_reporting import (
    CATEGORY_ORDER,
    G003_CATEGORIES,
    G004_CATEGORIES,
    G005_CATEGORIES,
    G003_SUMMARY_CSV,
    G004_SUMMARY_CSV,
    G005_SUMMARY_CSV,
    ROBUSTNESS_SUMMARY_CSV,
    FULL_DECISION_CSV,
    read_csv,
    bounded_source,
    select_current_bounded_row,
)
from build_scoring_models import _score_distribution_metrics
from run_g005_remaining_category_improvement import evaluation_negatives as g005_evaluation_negatives
from score_robustness_validation import AUDITED_HARD_NEGATIVE_REGIME, category_source_index, negative_regimes_for_score_target
from validate_subtyping_reason import build_entries, build_property_matched_negatives, load_category_smiles, sample_list

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
RESULTS_DIR = ROOT_DIR / 'results'
QED_RESULTS_DIR = RESULTS_DIR / 'qed_baseline_comparison'
QED_CSV = QED_RESULTS_DIR / 'qed_baseline_comparison.csv'
QED_JSON = QED_RESULTS_DIR / 'qed_baseline_comparison.json'
QED_SUMMARY_JSON = QED_RESULTS_DIR / 'qed_baseline_comparison_summary.json'

ROBUSTNESS_PER_RUN_CSV = RESULTS_DIR / 'score_robustness' / 'score_robustness_per_run.csv'
G003_PER_RUN_CSV = RESULTS_DIR / 'g003_broad_category_improvement' / 'g003_broad_category_improvement_per_run.csv'
G004_PER_RUN_CSV = RESULTS_DIR / 'g004_broad_drug_improvement' / 'g004_broad_category_improvement_per_run.csv'
G005_PER_RUN_CSV = RESULTS_DIR / 'g005_remaining_category_improvement' / 'g005_remaining_category_improvement_per_run.csv'


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: Any) -> str:
    if value in (None, ''):
        return ''
    return f"{float(value):.4f}".rstrip('0').rstrip('.')


def _per_run_row_by_source(category: str, current_row: dict[str, str] | None, source: str) -> dict[str, str] | None:
    if current_row is None:
        return None
    per_run_path = {
        'results/score_robustness/score_robustness_summary.csv': ROBUSTNESS_PER_RUN_CSV,
        'results/g003_broad_category_improvement/g003_broad_category_improvement_summary.csv': G003_PER_RUN_CSV,
        'results/g004_broad_drug_improvement/g004_broad_category_improvement_summary.csv': G004_PER_RUN_CSV,
        'results/g005_remaining_category_improvement/g005_remaining_category_improvement_summary.csv': G005_PER_RUN_CSV,
    }.get(source)
    if per_run_path is None or not per_run_path.exists():
        return None
    rows = read_csv(per_run_path)
    if source == 'results/score_robustness/score_robustness_summary.csv':
        regime = current_row.get('regime', '')
        for row in rows:
            if row.get('category') == category and row.get('regime') == regime:
                return row
        return None
    variant = current_row.get('variant', '')
    for row in rows:
        if row.get('category') == category and row.get('variant') == variant:
            return row
    return None


def _eval_negatives_for_current_row(
    *,
    slug: str,
    seed: int,
    positive_entries: list[dict[str, Any]],
    categories: dict[str, list[str]],
    category_entries: dict[str, list[dict[str, Any]]],
    all_other_entries: list[dict[str, Any]],
    source_index: dict[str, set[str]],
    eval_regime: str,
    unrepaired_negative_limit: int,
) -> list[str]:
    if eval_regime == AUDITED_HARD_NEGATIVE_REGIME:
        property_matched = build_property_matched_negatives(positive_entries, all_other_entries, unrepaired_negative_limit)
        negatives, _audit = g005_evaluation_negatives(
            slug=slug,
            positive_entries=positive_entries,
            all_target_smiles=set(categories[slug]),
            all_other_entries=all_other_entries,
            property_matched=property_matched,
            source_index=source_index,
            max_negative=unrepaired_negative_limit,
        )
        return negatives
    if eval_regime == 'property_matched':
        return build_property_matched_negatives(positive_entries, all_other_entries, unrepaired_negative_limit)
    regimes = negative_regimes_for_score_target(slug, categories, category_entries, unrepaired_negative_limit, seed)
    if eval_regime not in regimes:
        raise ValueError(f'Unsupported eval regime {eval_regime!r} for {slug!r}')
    return regimes[eval_regime]


def _qed_scores(entries: list[dict[str, Any]]) -> np.ndarray:
    return np.asarray([float(QED.qed(entry['mol'])) for entry in entries], dtype=float)


def build_qed_comparison_rows() -> list[dict[str, str]]:
    categories = load_category_smiles()
    category_entries = {slug: build_entries(smiles) for slug, smiles in categories.items()}
    source_index = category_source_index(category_entries)

    robustness = {(row['category'], row.get('regime', '')): row for row in read_csv(ROBUSTNESS_SUMMARY_CSV)}
    summaries = {
        'g003': read_csv(G003_SUMMARY_CSV),
        'g004': read_csv(G004_SUMMARY_CSV),
        'g005': read_csv(G005_SUMMARY_CSV),
    }
    full_by_category = {row['category']: row for row in read_csv(FULL_DECISION_CSV)}

    comparison_rows: list[dict[str, str]] = []
    for category in CATEGORY_ORDER:
        full = full_by_category.get(category)
        if full is None:
            continue
        current, _baseline, default_regime = select_current_bounded_row(category, summaries, robustness)
        source = bounded_source(category)
        per_run = _per_run_row_by_source(category, current, source)
        if per_run is None:
            comparison_rows.append(
                {
                    'output_id': category,
                    'output_name': full['name'],
                    'output_level': 'category',
                    'validation_regime': default_regime,
                    'sampled_positive_count': '',
                    'sampled_negative_count': '',
                    'current_reported_auc': fmt(current.get('auc_mean') if current else full.get('auc')),
                    'current_reported_balanced_accuracy': fmt(current.get('balanced_accuracy_mean') if current else full.get('balanced_accuracy')),
                    'qed_auc': '',
                    'qed_balanced_accuracy': '',
                    'delta_auc_vs_qed': '',
                    'delta_balanced_accuracy_vs_qed': '',
                    'beats_qed_auc': '',
                    'beats_qed_balanced_accuracy': '',
                    'reportable_vs_qed': '',
                    'comparison_note': 'No current per-run artifact was available to reconstruct the evaluation dataset for QED.',
                }
            )
            continue

        seed = int(per_run['seed'])
        positive_count = int(per_run['positive_count'])
        positive_smiles = sample_list(categories[category], positive_count, seed + 1000)
        positive_entries = build_entries(positive_smiles)
        all_other_entries = [entry for other_slug, entries in category_entries.items() if other_slug != category for entry in entries]

        if source == 'results/score_robustness/score_robustness_summary.csv':
            eval_regime = per_run.get('regime', default_regime)
            unrepaired_negative_limit = int(per_run.get('negative_count') or 0)
        else:
            eval_regime = per_run.get('eval_regime') or per_run.get('regime') or default_regime
            if eval_regime == AUDITED_HARD_NEGATIVE_REGIME and category in (G003_CATEGORIES | G004_CATEGORIES | G005_CATEGORIES):
                audit_path = {
                    'results/g003_broad_category_improvement/g003_broad_category_improvement_summary.csv': RESULTS_DIR / 'g003_broad_category_improvement' / 'g003_broad_category_improvement_audit.csv',
                    'results/g004_broad_drug_improvement/g004_broad_category_improvement_summary.csv': RESULTS_DIR / 'g004_broad_drug_improvement' / 'g004_broad_category_improvement_audit.csv',
                    'results/g005_remaining_category_improvement/g005_remaining_category_improvement_summary.csv': RESULTS_DIR / 'g005_remaining_category_improvement' / 'g005_remaining_category_improvement_audit.csv',
                }[source]
                audit_row = next(row for row in read_csv(audit_path) if row['category'] == category and int(row['seed']) == seed)
                unrepaired_negative_limit = int(audit_row['unrepaired_negative_count_before_audit'])
            else:
                unrepaired_negative_limit = int(per_run.get('negative_count') or 0)

        negative_smiles = _eval_negatives_for_current_row(
            slug=category,
            seed=seed,
            positive_entries=positive_entries,
            categories=categories,
            category_entries=category_entries,
            all_other_entries=all_other_entries,
            source_index=source_index,
            eval_regime=eval_regime,
            unrepaired_negative_limit=unrepaired_negative_limit,
        )
        negative_entries = build_entries(negative_smiles)
        metrics = _score_distribution_metrics(_qed_scores(positive_entries), _qed_scores(negative_entries))

        current_auc = float(current.get('auc_mean') or current.get('auc') or full.get('auc') or 0) if current else float(full.get('auc') or 0)
        current_ba = float(current.get('balanced_accuracy_mean') or current.get('balanced_accuracy') or full.get('balanced_accuracy') or 0) if current else float(full.get('balanced_accuracy') or 0)
        delta_auc = current_auc - float(metrics['auc'])
        delta_ba = current_ba - float(metrics['balanced_accuracy'])
        beats_auc = delta_auc > 0
        beats_ba = delta_ba > 0
        comparison_rows.append(
            {
                'output_id': category,
                'output_name': full['name'],
                'output_level': 'category',
                'validation_regime': eval_regime,
                'sampled_positive_count': str(len(positive_entries)),
                'sampled_negative_count': str(len(negative_entries)),
                'current_reported_auc': fmt(current_auc),
                'current_reported_balanced_accuracy': fmt(current_ba),
                'qed_auc': fmt(metrics['auc']),
                'qed_balanced_accuracy': fmt(metrics['balanced_accuracy']),
                'delta_auc_vs_qed': fmt(delta_auc),
                'delta_balanced_accuracy_vs_qed': fmt(delta_ba),
                'beats_qed_auc': 'yes' if beats_auc else 'no',
                'beats_qed_balanced_accuracy': 'yes' if beats_ba else 'no',
                'reportable_vs_qed': 'yes' if beats_auc else 'no',
                'comparison_note': 'Current reported score is compared against raw RDKit QED on the same reconstructed positive/evaluation-negative set.',
            }
        )
    return comparison_rows


def build_qed_summary(rows: list[dict[str, str]]) -> dict[str, Any]:
    valid_rows = [row for row in rows if row['qed_auc']]
    return {
        'comparison_basis': 'A category is marked reportable_vs_qed when the current reported AUC is higher than raw RDKit QED AUC on the same reconstructed evaluation set.',
        'categories_compared': len(valid_rows),
        'reportable_vs_qed': [row['output_id'] for row in valid_rows if row['reportable_vs_qed'] == 'yes'],
        'not_better_than_qed': [row['output_id'] for row in valid_rows if row['reportable_vs_qed'] == 'no'],
    }


def write_outputs(rows: list[dict[str, str]]) -> None:
    QED_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(QED_CSV, rows)
    QED_JSON.write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding='utf-8')
    QED_SUMMARY_JSON.write_text(json.dumps(build_qed_summary(rows), indent=2, ensure_ascii=False), encoding='utf-8')


def main() -> None:
    rows = build_qed_comparison_rows()
    write_outputs(rows)
    print(f'Wrote {len(rows)} QED comparison rows')


if __name__ == '__main__':
    main()
