from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import QED

from build_scoring_models import GEN, _choi_prepare_rows, _score_distribution_metrics
from validate_subtyping_reason import read_smiles

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
REBUILD_INPUTS_DIR = APP_DIR / 'output' / 'final_category_rebuild' / 'inputs'
MODELS_DIR = APP_DIR / 'output' / 'final_category_rebuild' / 'models'
RESULTS_DIR = ROOT_DIR / 'results' / 'final_category_rebuild'
SUMMARY_CSV = RESULTS_DIR / 'final_category_rebuild_summary.csv'
SUMMARY_JSON = RESULTS_DIR / 'final_category_rebuild_summary.json'
QED_CSV = RESULTS_DIR / 'final_category_rebuild_qed_comparison.csv'
QED_JSON = RESULTS_DIR / 'final_category_rebuild_qed_comparison.json'
QED_SUMMARY_JSON = RESULTS_DIR / 'final_category_rebuild_qed_comparison_summary.json'

FINAL_CATEGORIES = [
    'animal_drugs',
    'human_drugs',
    'cosmetics',
    'flavoring_agents',
    'food_additives',
    'food_contact_substances',
    'fragrances',
    'pesticides',
    'solvents',
    'surfactants',
]

QED_NOTE = 'Rebuilt broad-category scorer is compared against raw RDKit QED on the same retained-cross-category evaluation set after Tanimoto filtering.'


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def _fmt(value: float) -> float:
    return round(float(value), 4)


def _qed_scores(rows: list[dict[str, Any]]) -> np.ndarray:
    scores: list[float] = []
    for row in rows:
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol is not None:
            scores.append(float(QED.qed(mol)))
    return np.asarray(scores, dtype=float)


def _filter_negatives_after_tanimoto(positive_smiles: list[str], negative_source_smiles: list[str], threshold: float) -> tuple[list[str], list[str]]:
    pos_fps = []
    pos_valid = []
    for smiles in positive_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            pos_fps.append(GEN.GetFingerprint(mol))
            pos_valid.append(smiles)

    neg_fps = []
    neg_valid = []
    for smiles in negative_source_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            neg_fps.append(GEN.GetFingerprint(mol))
            neg_valid.append(smiles)

    if not pos_fps:
        return pos_valid, neg_valid

    max_sim = np.zeros(len(neg_fps), dtype=float)
    for pos_fp in pos_fps:
        max_sim = np.maximum(max_sim, np.asarray(DataStructs.BulkTanimotoSimilarity(pos_fp, neg_fps), dtype=float))
    filtered_negative = [smiles for smiles, similarity in zip(neg_valid, max_sim) if similarity < threshold]
    return pos_valid, filtered_negative


def build_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for category in FINAL_CATEGORIES:
        model = json.loads((MODELS_DIR / f'final_{category}.json').read_text(encoding='utf-8'))
        positive_smiles = read_smiles(REBUILD_INPUTS_DIR / f'{category}__positive.csv')
        negative_source_smiles = read_smiles(REBUILD_INPUTS_DIR / f'{category}__negative_source.csv')
        pos_valid, filtered_negative = _filter_negatives_after_tanimoto(positive_smiles, negative_source_smiles, threshold=float(model.get('tanimoto_threshold', 0.3)))
        pos_rows = _choi_prepare_rows(pos_valid)
        neg_rows = _choi_prepare_rows(filtered_negative)
        qed_metrics = _score_distribution_metrics(_qed_scores(pos_rows), _qed_scores(neg_rows))
        model_auc = float(model['metrics']['auc'])
        model_ba = float(model['metrics']['balanced_accuracy'])
        delta_auc = model_auc - float(qed_metrics['auc'])
        delta_ba = model_ba - float(qed_metrics['balanced_accuracy'])
        rows.append(
            {
                'category': category,
                'positive_count': len(pos_rows),
                'negative_count': len(neg_rows),
                'model_auc': _fmt(model_auc),
                'model_balanced_accuracy': _fmt(model_ba),
                'qed_auc': _fmt(qed_metrics['auc']),
                'qed_balanced_accuracy': _fmt(qed_metrics['balanced_accuracy']),
                'delta_auc_vs_qed': _fmt(delta_auc),
                'delta_balanced_accuracy_vs_qed': _fmt(delta_ba),
                'beats_qed_auc': 'yes' if delta_auc > 0 else 'no',
                'beats_qed_balanced_accuracy': 'yes' if delta_ba > 0 else 'no',
                'reportable_vs_qed': 'yes' if delta_auc > 0 and delta_ba > 0 else 'no',
                'comparison_note': QED_NOTE,
            }
        )
    return rows


def enrich_summary(qed_rows: list[dict[str, Any]]) -> None:
    qed_by_category = {row['category']: row for row in qed_rows}
    summary_rows = _read_csv(SUMMARY_CSV)
    enriched: list[dict[str, Any]] = []
    for row in summary_rows:
        merged: dict[str, Any] = dict(row)
        qed = qed_by_category.get(row['category'])
        if qed:
            merged.update(
                {
                    'qed_auc': qed['qed_auc'],
                    'qed_balanced_accuracy': qed['qed_balanced_accuracy'],
                    'delta_auc_vs_qed': qed['delta_auc_vs_qed'],
                    'delta_balanced_accuracy_vs_qed': qed['delta_balanced_accuracy_vs_qed'],
                    'beats_qed_auc': qed['beats_qed_auc'],
                    'beats_qed_balanced_accuracy': qed['beats_qed_balanced_accuracy'],
                    'reportable_vs_qed': qed['reportable_vs_qed'],
                }
            )
        enriched.append(merged)
    _write_csv(SUMMARY_CSV, enriched)
    SUMMARY_JSON.write_text(json.dumps(enriched, indent=2, ensure_ascii=False), encoding='utf-8')


def main() -> None:
    rows = build_rows()
    _write_csv(QED_CSV, rows)
    QED_JSON.write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding='utf-8')
    QED_SUMMARY_JSON.write_text(
        json.dumps(
            {
                'category_count': len(rows),
                'beats_qed_auc_count': sum(1 for row in rows if row['beats_qed_auc'] == 'yes'),
                'beats_qed_balanced_accuracy_count': sum(1 for row in rows if row['beats_qed_balanced_accuracy'] == 'yes'),
                'reportable_vs_qed_count': sum(1 for row in rows if row['reportable_vs_qed'] == 'yes'),
                'comparison_note': QED_NOTE,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding='utf-8',
    )
    enrich_summary(rows)
    print(f'Wrote {len(rows)} rebuilt-vs-QED comparison rows')


if __name__ == '__main__':
    main()
