from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
import optuna
from scipy.stats import ks_2samp
from sklearn.metrics import roc_auc_score, roc_curve

from build_scoring_models import (
    CHOI_CANDIDATE_PATTERNS,
    CHOI_DEFAULT_BAYES_TRIALS,
    CHOI_RANGE_OPTIONS,
    _choi_apply_category_priors,
    _choi_calc_props,

    _choi_precompute_pattern_rates,
    _choi_prepare_rows,
    _choi_step1_negative,
    _resolve_choi_category_prior,
    _write_model_config,
)
from validate_subtyping_reason import load_category_smiles, write_smiles_csv

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
OUTPUT_DIR = APP_DIR / 'output' / 'final_category_rebuild'
MODELS_DIR = OUTPUT_DIR / 'models'
INPUTS_DIR = OUTPUT_DIR / 'inputs'
RESULTS_DIR = ROOT_DIR / 'results' / 'final_category_rebuild'
SUMMARY_CSV = RESULTS_DIR / 'final_category_rebuild_summary.csv'
SUMMARY_JSON = RESULTS_DIR / 'final_category_rebuild_summary.json'

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

PROPERTY_KEYS = ['MW', 'logP', 'HBD', 'HBA', 'TPSA', 'RotBonds', 'FCsp3', 'AromaticRings']
GRID_KS = [0.05, 0.15, 0.25]
GRID_RATIO = [2.0, 3.0, 5.0]
GRID_W = [0.25, 0.4, 0.55, 0.7]

optuna.logging.set_verbosity(optuna.logging.WARNING)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    return value


def _collect_negative_source(categories: dict[str, list[str]], target: str) -> list[str]:
    positive = set(categories[target])
    negatives: list[str] = []
    seen: set[str] = set()
    for other in FINAL_CATEGORIES:
        if other == target:
            continue
        for smiles in categories[other]:
            if smiles in positive:
                continue
            if smiles in seen:
                continue
            seen.add(smiles)
            negatives.append(smiles)
    return negatives


def _rank_properties(pos_props, neg_props) -> list[tuple[str, float]]:
    ranks = []
    for feature in PROPERTY_KEYS:
        stat = float(ks_2samp(pos_props[feature].dropna(), neg_props[feature].dropna())[0])
        ranks.append((feature, stat))
    ranks.sort(key=lambda item: item[1], reverse=True)
    return ranks


def _select_props_with_fallback(pos_props, neg_props, ks_threshold: float, lower_q: float, upper_q: float) -> tuple[list[str], dict[str, tuple[float, float]]]:
    ranked = _rank_properties(pos_props, neg_props)
    selected = [name for name, stat in ranked if stat >= ks_threshold]
    if not selected:
        selected = [name for name, _ in ranked[: min(3, len(ranked))]]
    ranges = {
        feature: (
            float(pos_props[feature].quantile(lower_q)),
            float(pos_props[feature].quantile(upper_q)),
        )
        for feature in selected
    }
    return selected, ranges


def _select_patterns_with_fallback(pattern_rates: dict[str, dict[str, float]], ratio_threshold: float) -> dict[str, float]:
    selected: dict[str, float] = {}
    ranked: list[tuple[str, float, float, float]] = []
    for name, info in pattern_rates.items():
        pos_rate = float(info['pos_rate'])
        neg_rate = float(info['neg_rate'])
        ratio = pos_rate / neg_rate if neg_rate > 0 else (float('inf') if pos_rate > 0 else 0.0)
        delta = pos_rate - neg_rate
        ranked.append((name, ratio, delta, pos_rate))
        if ratio >= ratio_threshold and pos_rate > 0:
            selected[name] = delta
    if selected:
        return selected
    ranked.sort(key=lambda item: (item[1], item[2], item[3]), reverse=True)
    relaxed = [item for item in ranked if item[2] > 0 and item[3] > 0][:3]
    if not relaxed:
        relaxed = [item for item in ranked if item[3] > 0][:1]
    return {name: max(delta, 0.1) for name, _ratio, delta, _pos_rate in relaxed}
def _fast_distribution_metrics(pos_scores: np.ndarray, neg_scores: np.ndarray) -> dict[str, float]:
    if len(pos_scores) < 2 or len(neg_scores) < 2:
        return {
            'objective': 10.0,
            'auc': 0.0,
            'ks': 0.0,
            'balanced_accuracy': 0.0,
            'threshold': 0.5,
            'overlap': 1.0,
            'agro_hit_rate': 0.0,
            'negative_reject_rate': 0.0,
        }

    labels = np.r_[np.ones(len(pos_scores)), np.zeros(len(neg_scores))]
    scores = np.r_[pos_scores, neg_scores]
    try:
        auc = float(roc_auc_score(labels, scores))
    except ValueError:
        auc = 0.0

    ks_stat, _ = ks_2samp(pos_scores, neg_scores)
    fpr, tpr, thresholds = roc_curve(labels, scores)
    balanced_acc = 0.5 * (tpr + (1.0 - fpr))
    best_idx = int(np.argmax(balanced_acc))
    threshold_raw = float(thresholds[best_idx])
    threshold = threshold_raw if np.isfinite(threshold_raw) else 0.5
    threshold_penalty = abs(threshold - 0.5) if np.isfinite(threshold_raw) else 1.0
    best_balanced_acc = float(balanced_acc[best_idx])

    bins = np.linspace(-0.05, 1.05, 201)
    pos_density, _ = np.histogram(pos_scores, bins=bins, density=True)
    neg_density, _ = np.histogram(neg_scores, bins=bins, density=True)
    overlap_area = float(np.sum(np.minimum(pos_density, neg_density) * np.diff(bins)))

    objective = float(
        2.0 * overlap_area
        + (1.0 - auc)
        + (1.0 - best_balanced_acc)
        + 0.5 * (1.0 - float(ks_stat))
        + 0.2 * threshold_penalty
    )
    return {
        'objective': objective,
        'auc': auc,
        'ks': float(ks_stat),
        'balanced_accuracy': best_balanced_acc,
        'threshold': float(threshold),
        'overlap': overlap_area,
        'agro_hit_rate': float(np.mean(pos_scores >= threshold)),
        'negative_reject_rate': float(np.mean(neg_scores < threshold)),
    }


def _fast_metrics_for_w(pos_rows, neg_rows, selected_props: list[str], ranges: dict[str, tuple[float, float]], weights: dict[str, float], w: float) -> dict[str, float]:
    total_weight = sum(weights.values()) or 1.0

    def scores_for(rows) -> np.ndarray:
        scores = np.empty(len(rows), dtype=float)
        for idx, row in enumerate(rows):
            props = row['props']
            prop_score = sum(ranges[feature][0] <= props[feature] <= ranges[feature][1] for feature in selected_props) / len(selected_props)
            structure_score = sum(weight for name, weight in weights.items() if row['matches'][name]) / total_weight
            scores[idx] = w * prop_score + (1.0 - w) * structure_score
        return scores

    metrics = _fast_distribution_metrics(scores_for(pos_rows), scores_for(neg_rows))
    metrics['w'] = float(w)
    return metrics




def _grid_optimize(pos_rows, neg_rows, category_prior):
    pos_props = _choi_calc_props(pos_rows)
    neg_props = _choi_calc_props(neg_rows)
    pattern_rates = _choi_precompute_pattern_rates(pos_rows, neg_rows)
    best: dict[str, Any] | None = None
    for lower_q, upper_q in CHOI_RANGE_OPTIONS:
        for ks_threshold in GRID_KS:
            selected_props, ranges = _select_props_with_fallback(pos_props, neg_props, ks_threshold, lower_q, upper_q)
            for ratio_threshold in GRID_RATIO:
                pattern_weights = _select_patterns_with_fallback(pattern_rates, ratio_threshold)
                pattern_weights = _choi_apply_category_priors(pattern_weights, category_prior)
                if not selected_props or not pattern_weights:
                    continue
                for w in GRID_W:
                    metrics = _fast_metrics_for_w(pos_rows, neg_rows, selected_props, ranges, pattern_weights, w)
                    candidate = {
                        'selected_props': selected_props,
                        'ranges': ranges,
                        'pattern_weights': pattern_weights,
                        'metrics': metrics,
                        'lower_q': lower_q,
                        'upper_q': upper_q,
                        'ks_threshold': ks_threshold,
                        'ratio_threshold': ratio_threshold,
                        'optimization_method': 'mixed_grid_search',
                        'bayes_trials': 0,
                    }
                    if best is None or metrics['objective'] < best['metrics']['objective'] - 1e-12 or (
                        abs(metrics['objective'] - best['metrics']['objective']) <= 1e-12 and metrics['auc'] > best['metrics']['auc']
                    ):
                        best = candidate
    if best is None:
        raise ValueError('No mixed property+pattern candidate found.')
    return best


def _bayes_optimize(pos_rows, neg_rows, seed: int, bayes_trials: int, category_prior):
    pos_props = _choi_calc_props(pos_rows)
    neg_props = _choi_calc_props(neg_rows)
    pattern_rates = _choi_precompute_pattern_rates(pos_rows, neg_rows)

    def objective(trial: optuna.trial.Trial) -> float:
        range_idx = trial.suggest_categorical('range_idx', list(range(len(CHOI_RANGE_OPTIONS))))
        lower_q, upper_q = CHOI_RANGE_OPTIONS[range_idx]
        ks_threshold = trial.suggest_float('ks_threshold', 0.05, 0.35)
        ratio_threshold = trial.suggest_float('ratio_threshold', 1.5, 6.0)
        w = trial.suggest_float('w', 0.2, 0.8)
        selected_props, ranges = _select_props_with_fallback(pos_props, neg_props, ks_threshold, lower_q, upper_q)
        pattern_weights = _select_patterns_with_fallback(pattern_rates, ratio_threshold)
        pattern_weights = _choi_apply_category_priors(pattern_weights, category_prior)
        if not selected_props or not pattern_weights:
            return 10.0
        metrics = _fast_metrics_for_w(pos_rows, neg_rows, selected_props, ranges, pattern_weights, w)
        trial.set_user_attr('payload', {
            'selected_props': selected_props,
            'ranges': ranges,
            'pattern_weights': pattern_weights,
            'metrics': metrics,
            'lower_q': lower_q,
            'upper_q': upper_q,
            'ks_threshold': ks_threshold,
            'ratio_threshold': ratio_threshold,
        })
        return float(metrics['objective'])

    study = optuna.create_study(direction='minimize', sampler=optuna.samplers.TPESampler(seed=seed, multivariate=True))
    study.optimize(objective, n_trials=max(1, bayes_trials), show_progress_bar=False)
    payload = dict(study.best_trial.user_attrs.get('payload', {}))
    payload['optimization_method'] = 'mixed_lee_tpe_bayesian'
    payload['bayes_trials'] = int(max(1, bayes_trials))
    return payload


def rebuild_category(category: str, categories: dict[str, list[str]], seed: int, bayes_trials: int) -> tuple[dict[str, Any], dict[str, Any]]:
    positive_smiles = list(dict.fromkeys(categories[category]))
    negative_source_smiles = _collect_negative_source(categories, category)
    positive_csv = INPUTS_DIR / f'{category}__positive.csv'
    negative_source_csv = INPUTS_DIR / f'{category}__negative_source.csv'
    write_smiles_csv(positive_csv, positive_smiles)
    write_smiles_csv(negative_source_csv, negative_source_smiles)

    try:
        category_prior = _resolve_choi_category_prior('auto', category)
    except ValueError:
        category_prior = None

    pos_valid, filtered_negative = _choi_step1_negative(positive_smiles, negative_source_smiles, threshold=0.3)
    pos_rows = _choi_prepare_rows(pos_valid)
    neg_rows = _choi_prepare_rows(filtered_negative)
    if not pos_rows or not neg_rows:
        raise ValueError(f'Insufficient molecules after Tanimoto filtering for {category}.')

    grid = _grid_optimize(pos_rows, neg_rows, category_prior)
    bayes = _bayes_optimize(pos_rows, neg_rows, seed=seed, bayes_trials=bayes_trials, category_prior=category_prior)
    chosen = bayes if bayes.get('metrics', {}).get('objective', 10.0) <= grid['metrics']['objective'] + 1e-12 else grid

    config = {
        'model_id': f'final_{category}',
        'label': f"{category.replace('_', ' ').title()} / Final rebuilt scorer",
        'category': category,
        'model_type': 'choi_auto',
        'threshold': float(chosen['metrics'].get('threshold', 0.5)),
        'source_student': '최예빈 + 이서윤 hybrid final rebuild',
        'source_student_id': '20251266+20251288',
        'description': 'Final broad-category scorer rebuilt with retained cross-category negatives, positive-set physicochemical ranges, and structural patterns combined in one score.',
        'metrics': {
            'auc': float(chosen['metrics']['auc']),
            'ks': float(chosen['metrics']['ks']),
            'balanced_accuracy': float(chosen['metrics']['balanced_accuracy']),
            'overlap': float(chosen['metrics']['overlap']),
            'objective': float(chosen['metrics']['objective']),
            'negative_count': len(neg_rows),
        },
        'selected_props': chosen['selected_props'],
        'ranges': chosen['ranges'],
        'pattern_weights': chosen['pattern_weights'],
        'selected_patterns': {name: CHOI_CANDIDATE_PATTERNS[name] for name in chosen['pattern_weights']},
        'best_w': float(chosen['metrics'].get('w', 0.5)),
        'tanimoto_threshold': 0.3,
        'ks_threshold': float(chosen['ks_threshold']),
        'ratio_threshold': float(chosen['ratio_threshold']),
        'range_quantiles': [float(chosen['lower_q']), float(chosen['upper_q'])],
        'optimization_method': str(chosen['optimization_method']),
        'optimization_trials': int(chosen['bayes_trials']),
        'category_prior_used': category_prior is not None,
        'category_prior': category_prior,
        'positive_csv': str(positive_csv),
        'negative_source_csv': str(negative_source_csv),
        'negative_policy': {
            'final_category_pool': FINAL_CATEGORIES,
            'cross_category_overlap_policy': 'retain_overlap_from_other_categories',
            'excluded_only': [
                'exact_target_positive_overlap',
                'duplicate_smiles_in_negative_source',
                'near_positive_tanimoto_after_step1',
            ],
            'note': 'Molecules found in other categories are retained as negatives unless they are also in the target positive set.',
        },
    }
    model_path = MODELS_DIR / f'final_{category}.json'
    _write_model_config(config, model_path)
    summary = {
        'category': category,
        'positive_count': len(pos_rows),
        'negative_source_count': len(negative_source_smiles),
        'negative_count_after_tanimoto': len(neg_rows),
        'selected_prop_count': len(config['selected_props']),
        'selected_props': ';'.join(config['selected_props']),
        'pattern_count': len(config['selected_patterns']),
        'selected_patterns': ';'.join(config['selected_patterns'].keys()),
        'best_w': round(float(config['best_w']), 4),
        'auc': round(float(config['metrics']['auc']), 4),
        'ks': round(float(config['metrics']['ks']), 4),
        'balanced_accuracy': round(float(config['metrics']['balanced_accuracy']), 4),
        'threshold': round(float(config['threshold']), 4),
        'optimization_method': config['optimization_method'],
        'category_prior_used': config['category_prior_used'],
        'model_path': str(model_path),
    }
    return config, summary


def main() -> None:
    categories = load_category_smiles()
    missing = [category for category in FINAL_CATEGORIES if category not in categories]
    if missing:
        raise KeyError(f'Missing categories in input pool: {missing}')
    summaries: list[dict[str, Any]] = []
    configs: list[dict[str, Any]] = []
    for category in FINAL_CATEGORIES:
        config, summary = rebuild_category(category, categories, seed=42, bayes_trials=8)
        configs.append(config)
        summaries.append(summary)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(SUMMARY_CSV, summaries)
    SUMMARY_JSON.write_text(json.dumps(_json_ready(summaries), indent=2, ensure_ascii=False), encoding='utf-8')
    print(f'Rebuilt {len(summaries)} final category scorers')


if __name__ == '__main__':
    main()
