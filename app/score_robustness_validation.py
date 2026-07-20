from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
import statistics
from pathlib import Path
from typing import Any

from rdkit import DataStructs

from build_scoring_models import build_choi_model
from validate_subtyping_reason import (
    TMP_DIR,
    build_entries,
    build_property_matched_negatives,
    load_category_smiles,
    sample_list,
    write_smiles_csv,
)

APP_DIR = Path(__file__).resolve().parent
RESULTS_DIR = APP_DIR.parent / "results" / "score_robustness"
RUN_DIR = APP_DIR / "output" / "score_robustness"

SCORE_TARGET_SPECS: dict[str, dict[str, Any]] = {
    "pfas": {"related": ["surfactants", "solvents", "lipids", "food_contact_substances"]},
    "solvents": {"related": ["food_contact_substances", "food_additives", "flavoring_agents", "cosmetics"]},
    "flavoring_agents": {"related": ["fragrances", "food_additives", "cosmetics", "solvents"]},
    "food_additives": {"related": ["flavoring_agents", "food_contact_substances", "cosmetics", "human_drugs"]},
    "lipids": {"related": ["surfactants", "food_contact_substances", "cosmetics", "plasticizers"]},
    "endocrine_disruptors": {"related": ["human_drugs", "animal_drugs", "pesticides", "food_additives"]},
    "pesticides": {"related": ["endocrine_disruptors", "human_drugs", "animal_drugs", "food_additives"]},
    "fragrances": {"related": ["flavoring_agents", "cosmetics", "food_additives", "solvents"]},
    "surfactants": {"related": ["lipids", "cosmetics", "food_contact_substances", "solvents"]},
    "cosmetics": {"related": ["fragrances", "surfactants", "food_additives", "flavoring_agents", "solvents", "lipids", "plasticizers"]},
    "food_contact_substances": {"related": ["food_additives", "solvents", "lipids", "surfactants", "plasticizers", "cosmetics", "flavoring_agents"]},
    "human_drugs": {"related": ["animal_drugs", "endocrine_disruptors", "food_additives", "flavoring_agents", "cosmetics", "solvents", "pesticides"]},
    "animal_drugs": {"related": ["human_drugs", "endocrine_disruptors", "food_additives", "cosmetics", "solvents", "pesticides", "flavoring_agents"]},
}

EVIDENCE_PANEL_CATEGORIES = {"cosmetics", "food_contact_substances", "human_drugs", "animal_drugs"}
AUDITED_HARD_NEGATIVE_REGIME = "audited_hard_negative"
AUDITED_HARD_NEGATIVE_CATEGORIES = {
    "cosmetics",
    "food_contact_substances",
    "human_drugs",
    "animal_drugs",
    "flavoring_agents",
    "fragrances",
    "lipids",
}
NEAR_POSITIVE_TANIMOTO_THRESHOLD = 0.85
AUDITED_RETAINED_SET_TANIMOTO_THRESHOLD = 1.01
DUAL_USE_EXCLUSION_BY_TARGET = {
    "cosmetics": {"fragrances", "flavoring_agents", "surfactants", "lipids", "food_additives", "food_contact_substances"},
    "food_contact_substances": {"food_additives", "flavoring_agents", "lipids", "plasticizers", "cosmetics"},
    "human_drugs": {"animal_drugs", "endocrine_disruptors", "food_additives", "flavoring_agents", "cosmetics"},
    "animal_drugs": {"human_drugs", "endocrine_disruptors", "food_additives", "flavoring_agents", "cosmetics"},
    "flavoring_agents": {"fragrances", "food_additives", "cosmetics", "solvents"},
    "fragrances": {"flavoring_agents", "food_additives", "cosmetics", "solvents"},
    "lipids": {"surfactants", "food_contact_substances", "food_additives", "cosmetics", "plasticizers"},
}


def category_source_index(category_entries: dict[str, list[dict[str, Any]]]) -> dict[str, set[str]]:
    source_index: dict[str, set[str]] = {}
    for category, entries in category_entries.items():
        for entry in entries:
            source_index.setdefault(str(entry["smiles"]), set()).add(category)
    return source_index


def entry_index(entries: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for entry in entries:
        indexed.setdefault(str(entry["smiles"]), entry)
    return indexed


def build_audited_hard_negatives(
    slug: str,
    positive_entries: list[dict[str, Any]],
    all_target_smiles: set[str],
    property_matched_smiles: list[str],
    candidate_entry_index: dict[str, dict[str, Any]],
    source_index: dict[str, set[str]],
    max_negative: int,
) -> tuple[list[str], dict[str, Any]]:
    positive_smiles = {str(entry["smiles"]) for entry in positive_entries}
    positive_fps = [entry["fp"] for entry in positive_entries]
    dual_use_categories = DUAL_USE_EXCLUSION_BY_TARGET.get(slug, set())
    reason_counts: Counter[str] = Counter()
    retained: list[str] = []
    seen: set[str] = set()

    for smiles in property_matched_smiles:
        if smiles in seen:
            reason_counts["duplicate_candidate_smiles"] += 1
            continue
        seen.add(smiles)

        source_categories = source_index.get(smiles, set())
        if smiles in positive_smiles:
            reason_counts["exact_positive_overlap"] += 1
            continue
        if smiles in all_target_smiles:
            reason_counts["target_label_overlap"] += 1
            continue
        if source_categories & dual_use_categories:
            reason_counts["category_specific_dual_use"] += 1
            continue

        candidate_entry = candidate_entry_index.get(smiles)
        if candidate_entry is not None and positive_fps:
            best_similarity = max(DataStructs.BulkTanimotoSimilarity(candidate_entry["fp"], positive_fps))
            if best_similarity >= NEAR_POSITIVE_TANIMOTO_THRESHOLD:
                reason_counts["near_positive_tanimoto"] += 1
                continue

        retained.append(smiles)
        if len(retained) >= max_negative:
            break

    reason_counts["retained_audited_hard_negatives"] = len(retained)
    audit = {
        "category": slug,
        "regime": AUDITED_HARD_NEGATIVE_REGIME,
        "property_matched_pool_count": len(property_matched_smiles),
        "property_matched_unique_count": len(seen),
        "near_positive_tanimoto_threshold": NEAR_POSITIVE_TANIMOTO_THRESHOLD,
        "dual_use_categories": sorted(dual_use_categories),
        "reason_counts": dict(sorted(reason_counts.items())),
        "retained_count": len(retained),
    }
    return retained, audit



def select_positive_smiles(pool: list[str], limit: int, seed: int) -> list[str]:
    return sample_list(pool, limit, seed)


def negative_regimes_for_score_target(
    slug: str,
    categories: dict[str, list[str]],
    category_entries: dict[str, list[dict[str, Any]]],
    max_negative: int,
    seed: int,
) -> dict[str, list[str]]:
    related = set(SCORE_TARGET_SPECS[slug]["related"])

    all_other_smiles = []
    related_smiles = []
    distant_smiles = []
    candidate_entries = []

    for other_slug, smiles in categories.items():
        if other_slug == slug:
            continue
        all_other_smiles.extend(smiles)
        candidate_entries.extend(category_entries.get(other_slug, []))
        if other_slug in related:
            related_smiles.extend(smiles)
        else:
            distant_smiles.extend(smiles)

    all_other_smiles = list(dict.fromkeys(all_other_smiles))
    related_smiles = list(dict.fromkeys(related_smiles))
    distant_smiles = list(dict.fromkeys(distant_smiles))
    property_matched = build_property_matched_negatives(category_entries[slug], candidate_entries, max_negative)

    return {
        "all_other_random": sample_list(all_other_smiles, max_negative, seed + 11),
        "related_hard": sample_list(related_smiles, max_negative, seed + 22),
        "distant_other": sample_list(distant_smiles, max_negative, seed + 33),
        "property_matched": property_matched,
    }


def build_validation_model(
    slug: str,
    regime: str,
    seed: int,
    positive_smiles: list[str],
    negative_smiles: list[str],
    bayes_trials: int,
) -> dict[str, Any]:
    run_dir = RUN_DIR / slug / regime / f"seed_{seed}"
    positive_csv = run_dir / "positive.csv"
    negative_csv = run_dir / "negative.csv"
    model_path = run_dir / "model.json"
    write_smiles_csv(positive_csv, positive_smiles)
    write_smiles_csv(negative_csv, negative_smiles)
    tanimoto_threshold = (
        AUDITED_RETAINED_SET_TANIMOTO_THRESHOLD
        if regime == AUDITED_HARD_NEGATIVE_REGIME
        else 0.3
    )
    build_choi_model(
        positive_csv=positive_csv,
        negative_source_csv=negative_csv,
        model_id=f"robust_{slug}_{regime}_{seed}",
        label=f"{slug} {regime} seed {seed}",
        category=slug,
        output_path=model_path,
        tanimoto_threshold=tanimoto_threshold,
        bayes_trials=bayes_trials,
        seed=seed,
    )
    return json.loads(model_path.read_text(encoding="utf-8"))


def summarize_metrics(values: list[float]) -> dict[str, float]:
    return {
        "mean": round(statistics.mean(values), 4),
        "sd": round(statistics.stdev(values), 4) if len(values) > 1 else 0.0,
        "min": round(min(values), 4),
        "max": round(max(values), 4),
    }

def build_category_inventory(summary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in summary_rows:
        grouped.setdefault(str(row["category"]), []).append(row)

    inventory_rows = []
    for category, rows in grouped.items():
        rows_sorted = sorted(rows, key=lambda item: float(item["auc_mean"]))
        worst = rows_sorted[0]
        best = rows_sorted[-1]
        property_row = next((row for row in rows if row["regime"] == "property_matched"), worst)
        audited_row = next((row for row in rows if row["regime"] == AUDITED_HARD_NEGATIVE_REGIME), None)
        preferred_hard_negative_row = audited_row or property_row
        if category in EVIDENCE_PANEL_CATEGORIES:
            recommended_mode = "evidence_panel"
            retention_reason = "Keep in app as tracked experimental category with evidence-panel workflow despite weak or unstable global score."
        elif float(worst["auc_mean"]) >= 0.8 and float(worst["balanced_accuracy_mean"]) >= 0.75:
            recommended_mode = "robust_score"
            retention_reason = "Keep as release-grade scorer; performance remains strong across negative-set regimes."
        elif float(worst["auc_mean"]) >= 0.7:
            recommended_mode = "conditional_score"
            retention_reason = "Keep as score category, but mark robustness-sensitive and avoid overclaiming beyond validated regimes."
        else:
            recommended_mode = "experimental_score"
            retention_reason = "Keep in app and paper as experimental tracked category for future improvement, not as release-grade scorer."
        inventory_rows.append(
            {
                "category": category,
                "recommended_mode": recommended_mode,
                "kept_in_app": True,
                "worst_regime": worst["regime"],
                "worst_auc_mean": worst["auc_mean"],
                "worst_balanced_accuracy_mean": worst["balanced_accuracy_mean"],
                "best_regime": best["regime"],
                "best_auc_mean": best["auc_mean"],
                "property_matched_auc_mean": property_row["auc_mean"],
                "property_matched_balanced_accuracy_mean": property_row["balanced_accuracy_mean"],
                "preferred_hard_negative_regime": preferred_hard_negative_row["regime"],
                "preferred_hard_negative_auc_mean": preferred_hard_negative_row["auc_mean"],
                "preferred_hard_negative_balanced_accuracy_mean": preferred_hard_negative_row["balanced_accuracy_mean"],
                "retention_reason": retention_reason,
            }
        )
    inventory_rows.sort(key=lambda row: row["category"])
    return inventory_rows


def run_validation(
    seeds: list[int],
    max_positive: int,
    max_negative: int,
    bayes_trials: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    categories = load_category_smiles()
    category_entries = {slug: build_entries(smiles) for slug, smiles in categories.items()}
    source_index = category_source_index(category_entries)
    per_run_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    audit_run_rows: list[dict[str, Any]] = []
    audit_summary_rows: list[dict[str, Any]] = []

    for slug in SCORE_TARGET_SPECS:
        regime_names = ["all_other_random", "related_hard", "distant_other", "property_matched"]
        if slug in AUDITED_HARD_NEGATIVE_CATEGORIES:
            regime_names.append(AUDITED_HARD_NEGATIVE_REGIME)
        regime_metrics: dict[str, list[dict[str, float]]] = {regime: [] for regime in regime_names}
        for seed in seeds:
            positive_smiles = select_positive_smiles(categories[slug], max_positive, seed + 1000)
            positive_entries = build_entries(positive_smiles)
            regimes = negative_regimes_for_score_target(slug, categories, category_entries, max_negative=max_negative, seed=seed)
            original_entries = category_entries[slug]
            category_entries[slug] = positive_entries
            try:
                audit_row: dict[str, Any] | None = None
                seed_metric_rows: dict[str, dict[str, Any]] = {}
                all_other_entries = []
                for other_slug, entries in category_entries.items():
                    if other_slug != slug:
                        all_other_entries.extend(entries)
                property_matched = build_property_matched_negatives(positive_entries, all_other_entries, max_negative)
                regimes["property_matched"] = property_matched
                if slug in AUDITED_HARD_NEGATIVE_CATEGORIES:
                    audited_negatives, audit = build_audited_hard_negatives(
                        slug=slug,
                        positive_entries=positive_entries,
                        all_target_smiles=set(categories[slug]),
                        property_matched_smiles=property_matched,
                        candidate_entry_index=entry_index(all_other_entries),
                        source_index=source_index,
                        max_negative=max_negative,
                    )
                    regimes[AUDITED_HARD_NEGATIVE_REGIME] = audited_negatives
                    reason_counts = audit["reason_counts"]
                    audit_row = {
                        "category": slug,
                        "seed": seed,
                        "regime": AUDITED_HARD_NEGATIVE_REGIME,
                        "audit_schema_version": "g002_hard_negative_audit_v1",
                        "unrepaired_negative_count_before_audit": len(property_matched),
                        "repaired_negative_count_after_audit": len(audited_negatives),
                        "excluded_negative_count_by_audit": len(property_matched) - len(audited_negatives),
                        "property_matched_pool_count": audit["property_matched_pool_count"],
                        "property_matched_unique_count": audit["property_matched_unique_count"],
                        "retained_count": audit["retained_count"],
                        "duplicate_candidate_smiles": reason_counts.get("duplicate_candidate_smiles", 0),
                        "exact_positive_overlap": reason_counts.get("exact_positive_overlap", 0),
                        "target_label_overlap": reason_counts.get("target_label_overlap", 0),
                        "category_specific_dual_use": reason_counts.get("category_specific_dual_use", 0),
                        "near_positive_tanimoto": reason_counts.get("near_positive_tanimoto", 0),
                        "retained_audited_hard_negatives": reason_counts.get("retained_audited_hard_negatives", 0),
                        "post_audit_tanimoto_filter_applied": False,
                        "post_audit_tanimoto_threshold": "",
                        "near_positive_tanimoto_threshold": audit["near_positive_tanimoto_threshold"],
                        "dual_use_categories": ";".join(audit["dual_use_categories"]),
                        "unrepaired_metric_auc": "",
                        "unrepaired_metric_balanced_accuracy": "",
                        "unrepaired_metric_sd": "",
                        "repaired_metric_auc": "",
                        "repaired_metric_balanced_accuracy": "",
                        "repaired_metric_sd": "",
                        "delta_auc_repaired_minus_unrepaired": "",
                        "delta_ba_repaired_minus_unrepaired": "",
                    }
                    audit_run_rows.append(audit_row)

                for regime, negative_smiles in regimes.items():
                    config = build_validation_model(
                        slug=slug,
                        regime=regime,
                        seed=seed,
                        positive_smiles=positive_smiles,
                        negative_smiles=negative_smiles,
                        bayes_trials=bayes_trials,
                    )
                    metrics = config.get("metrics", {})
                    row = {
                        "category": slug,
                        "seed": seed,
                        "regime": regime,
                        "positive_count": len(positive_smiles),
                        "negative_count": int(metrics.get("negative_count", 0)),
                        "auc": round(float(metrics.get("auc", 0.0)), 4),
                        "ks": round(float(metrics.get("ks", 0.0)), 4),
                        "balanced_accuracy": round(float(metrics.get("balanced_accuracy", 0.0)), 4),
                        "overlap": round(float(metrics.get("overlap", 0.0)), 4),
                        "threshold": round(float(config.get("threshold", 0.5)), 4),
                        "optimization_method": str(config.get("optimization_method", "")),
                    }
                    per_run_rows.append(row)
                    regime_metrics[regime].append(row)
                    seed_metric_rows[regime] = row
                if audit_row is not None:
                    unrepaired_row = seed_metric_rows.get("property_matched")
                    repaired_row = seed_metric_rows.get(AUDITED_HARD_NEGATIVE_REGIME)
                    if unrepaired_row is not None and repaired_row is not None:
                        unrepaired_auc = float(unrepaired_row["auc"])
                        repaired_auc = float(repaired_row["auc"])
                        unrepaired_ba = float(unrepaired_row["balanced_accuracy"])
                        repaired_ba = float(repaired_row["balanced_accuracy"])
                        audit_row.update(
                            {
                                "unrepaired_metric_auc": unrepaired_row["auc"],
                                "unrepaired_metric_balanced_accuracy": unrepaired_row["balanced_accuracy"],
                                "repaired_metric_auc": repaired_row["auc"],
                                "repaired_metric_balanced_accuracy": repaired_row["balanced_accuracy"],
                                "delta_auc_repaired_minus_unrepaired": round(repaired_auc - unrepaired_auc, 4),
                                "delta_ba_repaired_minus_unrepaired": round(repaired_ba - unrepaired_ba, 4),
                            }
                        )
            finally:
                category_entries[slug] = original_entries

        for regime, rows in regime_metrics.items():
            auc_stats = summarize_metrics([float(row["auc"]) for row in rows])
            ks_stats = summarize_metrics([float(row["ks"]) for row in rows])
            ba_stats = summarize_metrics([float(row["balanced_accuracy"]) for row in rows])
            summary_rows.append(
                {
                    "category": slug,
                    "regime": regime,
                    "runs": len(rows),
                    "auc_mean": auc_stats["mean"],
                    "auc_sd": auc_stats["sd"],
                    "auc_min": auc_stats["min"],
                    "auc_max": auc_stats["max"],
                    "ks_mean": ks_stats["mean"],
                    "ks_sd": ks_stats["sd"],
                    "balanced_accuracy_mean": ba_stats["mean"],
                    "balanced_accuracy_sd": ba_stats["sd"],
                }
            )

    summary_by_category_regime = {(row["category"], row["regime"]): row for row in summary_rows}
    for category in sorted({row["category"] for row in audit_run_rows}):
        rows = [row for row in audit_run_rows if row["category"] == category]
        unrepaired_summary = summary_by_category_regime.get((category, "property_matched"), {})
        repaired_summary = summary_by_category_regime.get((category, AUDITED_HARD_NEGATIVE_REGIME), {})
        unrepaired_auc = float(unrepaired_summary.get("auc_mean", 0.0))
        repaired_auc = float(repaired_summary.get("auc_mean", 0.0))
        unrepaired_ba = float(unrepaired_summary.get("balanced_accuracy_mean", 0.0))
        repaired_ba = float(repaired_summary.get("balanced_accuracy_mean", 0.0))
        audit_summary_rows.append(
            {
                "category": category,
                "regime": AUDITED_HARD_NEGATIVE_REGIME,
                "runs": len(rows),
                "audit_schema_version": rows[0]["audit_schema_version"],
                "unrepaired_negative_count_before_audit": sum(int(row["unrepaired_negative_count_before_audit"]) for row in rows),
                "repaired_negative_count_after_audit": sum(int(row["repaired_negative_count_after_audit"]) for row in rows),
                "excluded_negative_count_by_audit": sum(int(row["excluded_negative_count_by_audit"]) for row in rows),
                "property_matched_pool_count": sum(int(row["property_matched_pool_count"]) for row in rows),
                "property_matched_unique_count": sum(int(row["property_matched_unique_count"]) for row in rows),
                "retained_count": sum(int(row["retained_count"]) for row in rows),
                "duplicate_candidate_smiles": sum(int(row["duplicate_candidate_smiles"]) for row in rows),
                "exact_positive_overlap": sum(int(row["exact_positive_overlap"]) for row in rows),
                "target_label_overlap": sum(int(row["target_label_overlap"]) for row in rows),
                "category_specific_dual_use": sum(int(row["category_specific_dual_use"]) for row in rows),
                "near_positive_tanimoto": sum(int(row["near_positive_tanimoto"]) for row in rows),
                "retained_audited_hard_negatives": sum(int(row["retained_audited_hard_negatives"]) for row in rows),
                "post_audit_tanimoto_filter_applied": False,
                "post_audit_tanimoto_threshold": "",
                "near_positive_tanimoto_threshold": NEAR_POSITIVE_TANIMOTO_THRESHOLD,
                "dual_use_categories": rows[0]["dual_use_categories"],
                "unrepaired_metric_auc": unrepaired_summary.get("auc_mean", ""),
                "unrepaired_metric_balanced_accuracy": unrepaired_summary.get("balanced_accuracy_mean", ""),
                "unrepaired_metric_sd": unrepaired_summary.get("auc_sd", ""),
                "repaired_metric_auc": repaired_summary.get("auc_mean", ""),
                "repaired_metric_balanced_accuracy": repaired_summary.get("balanced_accuracy_mean", ""),
                "repaired_metric_sd": repaired_summary.get("auc_sd", ""),
                "delta_auc_repaired_minus_unrepaired": round(repaired_auc - unrepaired_auc, 4),
                "delta_ba_repaired_minus_unrepaired": round(repaired_ba - unrepaired_ba, 4),
            }
        )
    return per_run_rows, summary_rows, audit_run_rows, audit_summary_rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run multi-seed robustness validation for score-based PubChem categories.")
    parser.add_argument("--seeds", type=int, nargs="+", default=[11, 23, 37])
    parser.add_argument("--max-positive", type=int, default=1000)
    parser.add_argument("--max-negative", type=int, default=3000)
    parser.add_argument("--bayes-trials", type=int, default=4)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    per_run_rows, summary_rows, audit_run_rows, audit_summary_rows = run_validation(
        seeds=args.seeds,
        max_positive=args.max_positive,
        max_negative=args.max_negative,
        bayes_trials=args.bayes_trials,
    )
    inventory_rows = build_category_inventory(summary_rows)

    write_csv(
        RESULTS_DIR / "score_robustness_per_run.csv",
        per_run_rows,
        [
            "category",
            "seed",
            "regime",
            "positive_count",
            "negative_count",
            "auc",
            "ks",
            "balanced_accuracy",
            "overlap",
            "threshold",
            "optimization_method",
        ],
    )
    write_csv(
        RESULTS_DIR / "score_robustness_summary.csv",
        summary_rows,
        [
            "category",
            "regime",
            "runs",
            "auc_mean",
            "auc_sd",
            "auc_min",
            "auc_max",
            "ks_mean",
            "ks_sd",
            "balanced_accuracy_mean",
            "balanced_accuracy_sd",
        ],
    )
    write_csv(
        RESULTS_DIR / "category_inventory.csv",
        inventory_rows,
        [
            "category",
            "recommended_mode",
            "kept_in_app",
            "worst_regime",
            "worst_auc_mean",
            "worst_balanced_accuracy_mean",
            "best_regime",
            "best_auc_mean",
            "property_matched_auc_mean",
            "property_matched_balanced_accuracy_mean",
            "preferred_hard_negative_regime",
            "preferred_hard_negative_auc_mean",
            "preferred_hard_negative_balanced_accuracy_mean",
            "retention_reason",
        ],
    )
    audit_fieldnames = [
        "category",
        "seed",
        "regime",
        "audit_schema_version",
        "unrepaired_negative_count_before_audit",
        "repaired_negative_count_after_audit",
        "excluded_negative_count_by_audit",
        "property_matched_pool_count",
        "property_matched_unique_count",
        "retained_count",
        "duplicate_candidate_smiles",
        "exact_positive_overlap",
        "target_label_overlap",
        "category_specific_dual_use",
        "near_positive_tanimoto",
        "retained_audited_hard_negatives",
        "post_audit_tanimoto_filter_applied",
        "post_audit_tanimoto_threshold",
        "near_positive_tanimoto_threshold",
        "dual_use_categories",
        "unrepaired_metric_auc",
        "unrepaired_metric_balanced_accuracy",
        "unrepaired_metric_sd",
        "repaired_metric_auc",
        "repaired_metric_balanced_accuracy",
        "repaired_metric_sd",
        "delta_auc_repaired_minus_unrepaired",
        "delta_ba_repaired_minus_unrepaired",
    ]
    write_csv(
        RESULTS_DIR / "audited_hard_negative_audit_per_run.csv",
        audit_run_rows,
        audit_fieldnames,
    )
    write_csv(
        RESULTS_DIR / "audited_hard_negative_audit_summary.csv",
        audit_summary_rows,
        ["category", "regime", "runs", *[field for field in audit_fieldnames if field not in {"category", "regime", "seed"}]],
    )
    (RESULTS_DIR / "score_robustness_summary.json").write_text(
        json.dumps(
            {
                "seeds": args.seeds,
                "summary": summary_rows,
                "inventory": inventory_rows,
                "audited_hard_negative_audit": audit_summary_rows,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(RESULTS_DIR)


if __name__ == "__main__":
    main()
