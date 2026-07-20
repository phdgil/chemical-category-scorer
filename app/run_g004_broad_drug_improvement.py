from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any

from build_scoring_models import (
    _choi_prepare_rows,
    _choi_two_scores,
    _score_distribution_metrics,
    build_choi_model,
)
from score_robustness_validation import (
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
RESULTS_DIR = ROOT_DIR / "results" / "g004_broad_drug_improvement"
RUN_DIR = APP_DIR / "output" / "g004_broad_drug_improvement"
SUBTYPE_RESULTS_DIR = ROOT_DIR / "results" / "official_drug_subtype_validation"
ARTIFACT_PREFIXES = ("g004_broad_drug_improvement", "g004_broad_category_improvement")

TARGET_CATEGORIES = ["human_drugs", "animal_drugs"]
VARIANTS: dict[str, dict[str, Any]] = {
    "baseline": {
        "train_regime": AUDITED_HARD_NEGATIVE_REGIME,
        "category_prior": None,
    },
    "category_prior": {
        "train_regime": AUDITED_HARD_NEGATIVE_REGIME,
        "category_prior": "auto",
    },
    "all_other_random": {
        "train_regime": "all_other_random",
        "category_prior": None,
    },
    "all_other_random_category_prior": {
        "train_regime": "all_other_random",
        "category_prior": "auto",
    },
    "related_hard": {
        "train_regime": "related_hard",
        "category_prior": None,
    },
    "related_hard_category_prior": {
        "train_regime": "related_hard",
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


def subtype_cue_sources() -> dict[str, Any]:
    sources = [
        SUBTYPE_RESULTS_DIR / "official_drug_subtype_partial_summary.json",
        SUBTYPE_RESULTS_DIR / "official_drug_subtype_animal_subset_summary.json",
        SUBTYPE_RESULTS_DIR / "official_drug_corrected_findings.json",
    ]
    return {
        "cue_policy": "official subtype outputs are recorded as coverage/error-analysis cues only; G004 trains and evaluates broad human_drugs/animal_drugs targets",
        "files_present": [str(path) for path in sources if path.exists()],
    }


def build_candidate_model(
    *,
    slug: str,
    variant: str,
    train_regime: str,
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
    tanimoto_threshold = (
        AUDITED_RETAINED_SET_TANIMOTO_THRESHOLD
        if train_regime == AUDITED_HARD_NEGATIVE_REGIME
        else 0.3
    )
    build_choi_model(
        positive_csv=positive_csv,
        negative_source_csv=negative_csv,
        model_id=f"g004_{slug}_{variant}_{seed}",
        label=f"G004 {slug} {variant} seed {seed}",
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
    config["evaluation_regime"] = "candidate_model_not_audited_evaluation"
    config["metric_regime"] = "training_regime_internal_model_metrics"
    config["artifact_semantics"] = {
        "story": "G004 broad human/animal drug improvement attempt",
        "artifact_role": "candidate_training_model",
        "target_scope": "broad_category_not_subtype",
        "training_negative_regime": train_regime,
        "training_negative_source_csv": str(negative_csv),
        "evaluation_regime": "audited_holdout_evaluation_is_external_to_model_json",
        "metric_regime": "model.metrics describe the training-regime model build, not audited evaluation metrics",
    }
    model_path.write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    return config, model_path


def evaluate_model_on_audited_holdout(
    config: dict[str, Any],
    positive_smiles: list[str],
    audited_negative_smiles: list[str],
) -> dict[str, float]:
    positive_rows = _choi_prepare_rows(positive_smiles)
    negative_rows = _choi_prepare_rows(audited_negative_smiles)
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
    positive_scores = best_w * pos_property_scores + (1.0 - best_w) * pos_structure_scores
    negative_scores = best_w * neg_property_scores + (1.0 - best_w) * neg_structure_scores
    metrics = _score_distribution_metrics(positive_scores, negative_scores)
    metrics["negative_count"] = len(negative_rows)
    return metrics


def training_negatives_for_variant(
    *,
    slug: str,
    variant: str,
    variant_spec: dict[str, Any],
    seed: int,
    categories: dict[str, list[str]],
    category_entries: dict[str, list[dict[str, Any]]],
    audited_negatives: list[str],
    max_negative: int,
) -> list[str]:
    train_regime = str(variant_spec["train_regime"])
    if train_regime == AUDITED_HARD_NEGATIVE_REGIME:
        return audited_negatives
    regimes = negative_regimes_for_score_target(
        slug,
        categories,
        category_entries,
        max_negative,
        seed + 4100,
    )
    if train_regime not in regimes:
        raise ValueError(f"Unsupported G004 train regime {train_regime!r} for {variant!r}.")
    return regimes[train_regime]


def run_experiment(
    *,
    seeds: list[int],
    max_positive: int,
    max_negative: int,
    bayes_trials: int,
) -> dict[str, Any]:
    categories = load_category_smiles()
    category_entries = {slug: build_entries(smiles) for slug, smiles in categories.items()}
    source_index = category_source_index(category_entries)

    per_run_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []

    for slug in TARGET_CATEGORIES:
        for seed in seeds:
            positive_smiles = sample_list(categories[slug], max_positive, seed + 1000)
            positive_entries = build_entries(positive_smiles)
            all_other_entries = []
            for other_slug, entries in category_entries.items():
                if other_slug != slug:
                    all_other_entries.extend(entries)
            property_matched = build_property_matched_negatives(positive_entries, all_other_entries, max_negative)
            audited_negatives, audit = build_audited_hard_negatives(
                slug=slug,
                positive_entries=positive_entries,
                all_target_smiles=set(categories[slug]),
                property_matched_smiles=property_matched,
                candidate_entry_index=entry_index(all_other_entries),
                source_index=source_index,
                max_negative=max_negative,
            )
            reason_counts = audit["reason_counts"]
            audit_rows.append(
                {
                    "category": slug,
                    "seed": seed,
                    "regime": AUDITED_HARD_NEGATIVE_REGIME,
                    "unrepaired_negative_count_before_audit": len(property_matched),
                    "repaired_negative_count_after_audit": len(audited_negatives),
                    "excluded_negative_count_by_audit": len(property_matched) - len(audited_negatives),
                    "property_matched_pool_count": audit["property_matched_pool_count"],
                    "property_matched_unique_count": audit["property_matched_unique_count"],
                    "retained_count": audit["retained_count"],
                    "category_specific_dual_use": reason_counts.get("category_specific_dual_use", 0),
                    "near_positive_tanimoto": reason_counts.get("near_positive_tanimoto", 0),
                    "retained_audited_hard_negatives": reason_counts.get("retained_audited_hard_negatives", 0),
                    "near_positive_tanimoto_threshold": audit["near_positive_tanimoto_threshold"],
                    "dual_use_categories": ";".join(audit["dual_use_categories"]),
                }
            )
            for variant, variant_spec in VARIANTS.items():
                train_regime = str(variant_spec["train_regime"])
                train_negatives = training_negatives_for_variant(
                    slug=slug,
                    variant=variant,
                    variant_spec=variant_spec,
                    seed=seed,
                    categories=categories,
                    category_entries=category_entries,
                    audited_negatives=audited_negatives,
                    max_negative=max_negative,
                )
                config, model_path = build_candidate_model(
                    slug=slug,
                    variant=variant,
                    train_regime=train_regime,
                    category_prior=variant_spec.get("category_prior"),
                    seed=seed,
                    positive_smiles=positive_smiles,
                    negative_smiles=train_negatives,
                    bayes_trials=bayes_trials,
                )
                metrics = evaluate_model_on_audited_holdout(config, positive_smiles, audited_negatives)
                category_prior = config.get("category_prior") or {}
                per_run_rows.append(
                    {
                        "category": slug,
                        "variant": variant,
                        "seed": seed,
                        "regime": AUDITED_HARD_NEGATIVE_REGIME,
                        "train_regime": train_regime,
                        "eval_regime": AUDITED_HARD_NEGATIVE_REGIME,
                        "target_scope": "broad_category",
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
        baseline_rows = [row for row in per_run_rows if row["category"] == slug and row["variant"] == "baseline"]
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
                    "regime": AUDITED_HARD_NEGATIVE_REGIME,
                    "train_regime": rows[0]["train_regime"],
                    "eval_regime": AUDITED_HARD_NEGATIVE_REGIME,
                    "target_scope": "broad_category",
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
            "train_eval_separation": "per-run audited metrics use eval_regime; model JSON metrics remain training-regime build metrics",
        },
        "subtype_cue_sources": subtype_cue_sources(),
        "per_run": per_run_rows,
        "summary": summary_rows,
        "audit": audit_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run bounded G004 broad human/animal drug improvement attempts.")
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
        "regime",
        "train_regime",
        "eval_regime",
        "target_scope",
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
        "regime",
        "train_regime",
        "eval_regime",
        "target_scope",
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
        "regime",
        "unrepaired_negative_count_before_audit",
        "repaired_negative_count_after_audit",
        "excluded_negative_count_by_audit",
        "property_matched_pool_count",
        "property_matched_unique_count",
        "retained_count",
        "category_specific_dual_use",
        "near_positive_tanimoto",
        "retained_audited_hard_negatives",
        "near_positive_tanimoto_threshold",
        "dual_use_categories",
    ]
    for artifact_prefix in ARTIFACT_PREFIXES:
        write_csv(
            RESULTS_DIR / f"{artifact_prefix}_per_run.csv",
            payload["per_run"],
            per_run_fieldnames,
        )
        write_csv(
            RESULTS_DIR / f"{artifact_prefix}_summary.csv",
            payload["summary"],
            summary_fieldnames,
        )
        write_csv(
            RESULTS_DIR / f"{artifact_prefix}_audit.csv",
            payload["audit"],
            audit_fieldnames,
        )
        (RESULTS_DIR / f"{artifact_prefix}.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    print(RESULTS_DIR)


if __name__ == "__main__":
    main()
