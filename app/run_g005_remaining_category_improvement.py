from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any

from build_scoring_models import (
    CHOI_CATEGORY_PRIORS,
    LEE_DIR,
    _choi_prepare_rows,
    _choi_two_scores,
    _score_distribution_metrics,
    build_choi_model,
)
from score_robustness_validation import (
    AUDITED_HARD_NEGATIVE_CATEGORIES,
    AUDITED_HARD_NEGATIVE_REGIME,
    AUDITED_RETAINED_SET_TANIMOTO_THRESHOLD,
    build_audited_hard_negatives,
    category_source_index,
    entry_index,
    negative_regimes_for_score_target,
)
from validate_subtyping_reason import (
    TMP_DIR,
    build_entries,
    build_property_matched_negatives,
    load_category_smiles,
    sample_list,
    write_smiles_csv,
)

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
RESULTS_DIR = ROOT_DIR / "results" / "g005_remaining_category_improvement"
RUN_DIR = APP_DIR / "output" / "g005_remaining_category_improvement"

TARGET_CATEGORIES = [
    "pesticides",
    "fragrances",
    "food_additives",
    "flavoring_agents",
    "lipids",
    "solvents",
    "endocrine_disruptors",
    "surfactants",
]

PRIOR_SUPPORTED_CATEGORIES = {
    slug for slug, spec in CHOI_CATEGORY_PRIORS.items() if "alias_for" not in spec
}

VARIANTS: dict[str, dict[str, Any]] = {
    "baseline_hard_negative": {"train_regime": "eval_hard_negative", "category_prior": None},
    "category_prior_hard_negative": {"train_regime": "eval_hard_negative", "category_prior": "auto", "requires_prior": True},
    "all_other_random": {"train_regime": "all_other_random", "category_prior": None},
    "all_other_random_category_prior": {"train_regime": "all_other_random", "category_prior": "auto", "requires_prior": True},
    "related_hard": {"train_regime": "related_hard", "category_prior": None},
    "related_hard_category_prior": {"train_regime": "related_hard", "category_prior": "auto", "requires_prior": True},
    "pesticides_native_lee_negatives": {
        "targets": {"pesticides"},
        "train_regime": "pesticides_native_lee_negative_data",
        "category_prior": None,
    },
    "fragrances_native_pesticide_negatives": {
        "targets": {"fragrances"},
        "train_regime": "choi_native_pesticides_negative_source",
        "category_prior": None,
    },
    "fragrances_native_pesticide_negatives_category_prior": {
        "targets": {"fragrances"},
        "train_regime": "choi_native_pesticides_negative_source",
        "category_prior": "auto",
    },
    "surfactants_native_pesticide_negatives": {
        "targets": {"surfactants"},
        "train_regime": "choi_native_pesticides_negative_source",
        "category_prior": None,
    },
    "surfactants_native_pesticide_negatives_category_prior": {
        "targets": {"surfactants"},
        "train_regime": "choi_native_pesticides_negative_source",
        "category_prior": "auto",
    },
}


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize(values: list[float]) -> dict[str, float]:
    return {
        "mean": round(statistics.mean(values), 4),
        "sd": round(statistics.stdev(values), 4) if len(values) > 1 else 0.0,
        "min": round(min(values), 4),
        "max": round(max(values), 4),
    }


def read_smiles_flexible(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        columns = ["SMILES", "Smiles", "smiles", "CanonicalSMILES", "canonical_smiles"]
        column = next((name for name in columns if name in reader.fieldnames), None)
        if column is None:
            raise KeyError(f"No SMILES column found in {path}")
        return list(dict.fromkeys(row[column].strip() for row in reader if row.get(column, "").strip()))


def build_candidate_model(
    *,
    slug: str,
    variant: str,
    train_regime: str,
    eval_regime: str,
    category_prior: str | None,
    seed: int,
    positive_smiles: list[str],
    negative_smiles: list[str],
    bayes_trials: int,
) -> tuple[dict[str, Any], Path]:
    run_dir = RUN_DIR / slug / variant / f"seed_{seed}"
    positive_csv = run_dir / "positive.csv"
    negative_csv = run_dir / f"{train_regime}.csv"
    model_path = run_dir / "model.json"
    write_smiles_csv(positive_csv, positive_smiles)
    write_smiles_csv(negative_csv, negative_smiles)
    tanimoto_threshold = AUDITED_RETAINED_SET_TANIMOTO_THRESHOLD if train_regime == AUDITED_HARD_NEGATIVE_REGIME else 0.3
    build_choi_model(
        positive_csv=positive_csv,
        negative_source_csv=negative_csv,
        model_id=f"g005_{slug}_{variant}_{seed}",
        label=f"G005 {slug} {variant} seed {seed}",
        category=slug,
        output_path=model_path,
        tanimoto_threshold=tanimoto_threshold,
        bayes_trials=bayes_trials,
        seed=seed,
        category_prior=category_prior,
    )
    config = json.loads(model_path.read_text(encoding="utf-8"))
    config["training_regime"] = train_regime
    config["training_negative_source"] = str(negative_csv)
    config["evaluation_regime"] = "candidate_model_not_external_evaluation"
    config["metric_regime"] = "training_regime_internal_model_metrics"
    config["artifact_semantics"] = {
        "story": "G005 remaining broad-category improvement attempt",
        "artifact_role": "candidate_training_model",
        "target_scope": "broad_category_not_subtype",
        "training_negative_regime": train_regime,
        "training_negative_source_csv": str(negative_csv),
        "external_evaluation_regime": eval_regime,
        "metric_regime": "model.metrics describe the training-regime model build; per-run CSV/JSON rows contain external evaluation metrics",
    }
    model_path.write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    return config, model_path


def evaluate_model(config: dict[str, Any], positive_smiles: list[str], negative_smiles: list[str]) -> dict[str, float]:
    positive_rows = _choi_prepare_rows(positive_smiles)
    negative_rows = _choi_prepare_rows(negative_smiles)
    pos_property_scores, pos_structure_scores = _choi_two_scores(
        positive_rows,
        config.get("selected_props", []),
        config.get("ranges", {}),
        config.get("pattern_weights", {}),
    )
    neg_property_scores, neg_structure_scores = _choi_two_scores(
        negative_rows,
        config.get("selected_props", []),
        config.get("ranges", {}),
        config.get("pattern_weights", {}),
    )
    best_w = float(config.get("best_w", 0.5))
    metrics = _score_distribution_metrics(
        best_w * pos_property_scores + (1.0 - best_w) * pos_structure_scores,
        best_w * neg_property_scores + (1.0 - best_w) * neg_structure_scores,
    )
    metrics["negative_count"] = len(negative_rows)
    return metrics


def evaluation_negatives(
    *,
    slug: str,
    positive_entries: list[dict[str, Any]],
    all_target_smiles: set[str],
    all_other_entries: list[dict[str, Any]],
    property_matched: list[str],
    source_index: dict[str, set[str]],
    max_negative: int,
) -> tuple[list[str], dict[str, Any]]:
    if slug in AUDITED_HARD_NEGATIVE_CATEGORIES:
        negatives, audit = build_audited_hard_negatives(
            slug=slug,
            positive_entries=positive_entries,
            all_target_smiles=all_target_smiles,
            property_matched_smiles=property_matched,
            candidate_entry_index=entry_index(all_other_entries),
            source_index=source_index,
            max_negative=max_negative,
        )
        audit["eval_regime"] = AUDITED_HARD_NEGATIVE_REGIME
        return negatives, audit
    return property_matched, {
        "eval_regime": "property_matched",
        "property_matched_pool_count": len(property_matched),
        "property_matched_unique_count": len(set(property_matched)),
        "retained_count": len(property_matched),
        "near_positive_tanimoto_threshold": "not_audited",
        "dual_use_categories": [],
        "reason_counts": {"retained_property_matched_negatives": len(property_matched)},
    }


def training_negatives_for_variant(
    *,
    slug: str,
    variant_spec: dict[str, Any],
    seed: int,
    categories: dict[str, list[str]],
    category_entries: dict[str, list[dict[str, Any]]],
    eval_negatives: list[str],
    max_negative: int,
) -> tuple[str, list[str]]:
    train_regime = str(variant_spec["train_regime"])
    if train_regime == "eval_hard_negative":
        return (AUDITED_HARD_NEGATIVE_REGIME if slug in AUDITED_HARD_NEGATIVE_CATEGORIES else "property_matched"), eval_negatives
    if train_regime == "choi_native_pesticides_negative_source":
        return train_regime, sample_list(categories["pesticides"], max_negative, seed + 3100)
    if train_regime == "pesticides_native_lee_negative_data":
        return train_regime, sample_list(read_smiles_flexible(LEE_DIR / "negative_data.csv"), max_negative, seed + 3200)
    regimes = negative_regimes_for_score_target(slug, categories, category_entries, max_negative, seed + 4100)
    if train_regime not in regimes:
        raise ValueError(f"Unsupported G005 train regime {train_regime!r} for {slug!r}.")
    return train_regime, regimes[train_regime]


def run_experiment(*, seeds: list[int], max_positive: int, max_negative: int, bayes_trials: int) -> dict[str, Any]:
    categories = load_category_smiles()
    category_entries = {slug: build_entries(smiles) for slug, smiles in categories.items()}
    source_index = category_source_index(category_entries)
    per_run_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []

    for slug in TARGET_CATEGORIES:
        for seed in seeds:
            positive_smiles = sample_list(categories[slug], max_positive, seed + 1000)
            positive_entries = build_entries(positive_smiles)
            all_other_entries = [entry for other_slug, entries in category_entries.items() if other_slug != slug for entry in entries]
            property_matched = build_property_matched_negatives(positive_entries, all_other_entries, max_negative)
            eval_negatives, audit = evaluation_negatives(
                slug=slug,
                positive_entries=positive_entries,
                all_target_smiles=set(categories[slug]),
                all_other_entries=all_other_entries,
                property_matched=property_matched,
                source_index=source_index,
                max_negative=max_negative,
            )
            reason_counts = audit["reason_counts"]
            eval_regime = str(audit["eval_regime"])
            audit_rows.append(
                {
                    "category": slug,
                    "seed": seed,
                    "eval_regime": eval_regime,
                    "unrepaired_negative_count_before_audit": len(property_matched),
                    "repaired_negative_count_after_audit": len(eval_negatives),
                    "excluded_negative_count_by_audit": len(property_matched) - len(eval_negatives),
                    "property_matched_pool_count": audit["property_matched_pool_count"],
                    "property_matched_unique_count": audit["property_matched_unique_count"],
                    "retained_count": audit["retained_count"],
                    "category_specific_dual_use": reason_counts.get("category_specific_dual_use", 0),
                    "near_positive_tanimoto": reason_counts.get("near_positive_tanimoto", 0),
                    "retained_eval_negatives": reason_counts.get("retained_audited_hard_negatives", reason_counts.get("retained_property_matched_negatives", 0)),
                    "near_positive_tanimoto_threshold": audit["near_positive_tanimoto_threshold"],
                    "dual_use_categories": ";".join(audit["dual_use_categories"]),
                }
            )
            for variant, variant_spec in VARIANTS.items():
                targets = variant_spec.get("targets")
                if targets is not None and slug not in targets:
                    continue
                if variant_spec.get("requires_prior") and slug not in PRIOR_SUPPORTED_CATEGORIES:
                    continue
                train_regime, train_negatives = training_negatives_for_variant(
                    slug=slug,
                    variant_spec=variant_spec,
                    seed=seed,
                    categories=categories,
                    category_entries=category_entries,
                    eval_negatives=eval_negatives,
                    max_negative=max_negative,
                )
                config, model_path = build_candidate_model(
                    slug=slug,
                    variant=variant,
                    train_regime=train_regime,
                    eval_regime=eval_regime,
                    category_prior=variant_spec.get("category_prior"),
                    seed=seed,
                    positive_smiles=positive_smiles,
                    negative_smiles=train_negatives,
                    bayes_trials=bayes_trials,
                )
                metrics = evaluate_model(config, positive_smiles, eval_negatives)
                category_prior = config.get("category_prior") or {}
                per_run_rows.append(
                    {
                        "category": slug,
                        "variant": variant,
                        "seed": seed,
                        "train_regime": train_regime,
                        "eval_regime": eval_regime,
                        "target_scope": "broad_category",
                        "metric_semantics": "external_eval_metrics_against_eval_regime_negatives",
                        "positive_count": len(positive_smiles),
                        "train_negative_source_count": len(train_negatives),
                        "negative_count": int(metrics.get("negative_count", 0)),
                        "auc": round(float(metrics.get("auc", 0.0)), 4),
                        "ks": round(float(metrics.get("ks", 0.0)), 4),
                        "balanced_accuracy": round(float(metrics.get("balanced_accuracy", 0.0)), 4),
                        "overlap": round(float(metrics.get("overlap", 0.0)), 4),
                        "threshold": round(float(metrics.get("threshold", config.get("threshold", 0.5))), 4),
                        "best_w": round(float(config.get("best_w", 0.5)), 4),
                        "optimization_method": str(config.get("optimization_method", "")),
                        "category_prior_used": bool(config.get("category_prior_used", False)),
                        "prior_source": str(category_prior.get("source", "")),
                        "prior_motifs": ";".join(category_prior.get("motifs", {}).keys()),
                        "model_path": str(model_path),
                    }
                )

    summary_rows: list[dict[str, Any]] = []
    for slug in TARGET_CATEGORIES:
        baseline_rows = [row for row in per_run_rows if row["category"] == slug and row["variant"] == "baseline_hard_negative"]
        baseline_auc_mean = statistics.mean(float(row["auc"]) for row in baseline_rows)
        baseline_ba_mean = statistics.mean(float(row["balanced_accuracy"]) for row in baseline_rows)
        for variant in VARIANTS:
            rows = [row for row in per_run_rows if row["category"] == slug and row["variant"] == variant]
            if not rows:
                continue
            auc_stats = summarize([float(row["auc"]) for row in rows])
            ba_stats = summarize([float(row["balanced_accuracy"]) for row in rows])
            ks_stats = summarize([float(row["ks"]) for row in rows])
            summary_rows.append(
                {
                    "category": slug,
                    "variant": variant,
                    "train_regime": rows[0]["train_regime"],
                    "eval_regime": rows[0]["eval_regime"],
                    "target_scope": "broad_category",
                    "metric_semantics": "external_eval_metrics_against_eval_regime_negatives",
                    "runs": len(rows),
                    "auc_mean": auc_stats["mean"],
                    "auc_sd": auc_stats["sd"],
                    "auc_min": auc_stats["min"],
                    "auc_max": auc_stats["max"],
                    "balanced_accuracy_mean": ba_stats["mean"],
                    "balanced_accuracy_sd": ba_stats["sd"],
                    "ks_mean": ks_stats["mean"],
                    "ks_sd": ks_stats["sd"],
                    "delta_auc_vs_baseline_mean": round(auc_stats["mean"] - baseline_auc_mean, 4),
                    "delta_ba_vs_baseline_mean": round(ba_stats["mean"] - baseline_ba_mean, 4),
                }
            )
    return {
        "parameters": {
            "seeds": seeds,
            "max_positive": max_positive,
            "max_negative": max_negative,
            "bayes_trials": bayes_trials,
            "targets": TARGET_CATEGORIES,
            "variants": list(VARIANTS),
            "target_scope": "broad_category",
            "eval_policy": "audited hard negatives for locally audited categories; property_matched for remaining broad categories",
        },
        "per_run": per_run_rows,
        "summary": summary_rows,
        "audit": audit_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run bounded G005 remaining broad-category improvement attempts.")
    parser.add_argument("--seeds", type=int, nargs="+", default=[11, 23, 37])
    parser.add_argument("--max-positive", type=int, default=1000)
    parser.add_argument("--max-negative", type=int, default=3000)
    parser.add_argument("--bayes-trials", type=int, default=4)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    payload = run_experiment(
        seeds=args.seeds,
        max_positive=args.max_positive,
        max_negative=args.max_negative,
        bayes_trials=args.bayes_trials,
    )
    per_run_fieldnames = [
        "category",
        "variant",
        "seed",
        "train_regime",
        "eval_regime",
        "target_scope",
        "metric_semantics",
        "positive_count",
        "train_negative_source_count",
        "negative_count",
        "auc",
        "ks",
        "balanced_accuracy",
        "overlap",
        "threshold",
        "best_w",
        "optimization_method",
        "category_prior_used",
        "prior_source",
        "prior_motifs",
        "model_path",
    ]
    summary_fieldnames = [
        "category",
        "variant",
        "train_regime",
        "eval_regime",
        "target_scope",
        "metric_semantics",
        "runs",
        "auc_mean",
        "auc_sd",
        "auc_min",
        "auc_max",
        "balanced_accuracy_mean",
        "balanced_accuracy_sd",
        "ks_mean",
        "ks_sd",
        "delta_auc_vs_baseline_mean",
        "delta_ba_vs_baseline_mean",
    ]
    audit_fieldnames = [
        "category",
        "seed",
        "eval_regime",
        "unrepaired_negative_count_before_audit",
        "repaired_negative_count_after_audit",
        "excluded_negative_count_by_audit",
        "property_matched_pool_count",
        "property_matched_unique_count",
        "retained_count",
        "category_specific_dual_use",
        "near_positive_tanimoto",
        "retained_eval_negatives",
        "near_positive_tanimoto_threshold",
        "dual_use_categories",
    ]
    write_csv(RESULTS_DIR / "g005_remaining_category_improvement_per_run.csv", payload["per_run"], per_run_fieldnames)
    write_csv(RESULTS_DIR / "g005_remaining_category_improvement_summary.csv", payload["summary"], summary_fieldnames)
    write_csv(RESULTS_DIR / "g005_remaining_category_improvement_audit.csv", payload["audit"], audit_fieldnames)
    (RESULTS_DIR / "g005_remaining_category_improvement.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(RESULTS_DIR)


if __name__ == "__main__":
    main()
