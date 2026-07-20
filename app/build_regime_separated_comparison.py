from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
RESULTS_DIR = ROOT_DIR / "results"
PAPER_DIR = ROOT_DIR / "paper"
DATA_MODELS_DIR = APP_DIR / "data" / "models"
OUTPUT_DIR = APP_DIR / "output"

OUT_CSV = RESULTS_DIR / "regime_separated_comparison_table.csv"
OUT_JSON = RESULTS_DIR / "regime_separated_comparison_table.json"
OUT_BASELINE_MATRIX_CSV = RESULTS_DIR / "score_robustness" / "baseline_evidence_matrix.csv"
OUT_BASELINE_MATRIX_JSON = RESULTS_DIR / "score_robustness" / "baseline_evidence_matrix.json"
OUT_PAPER_CSV = PAPER_DIR / "regime_separated_comparison_table.csv"

SCORE_SUMMARY_CSV = RESULTS_DIR / "score_robustness" / "score_robustness_summary.csv"
SCORE_PER_RUN_CSV = RESULTS_DIR / "score_robustness" / "score_robustness_per_run.csv"
AUDITED_HARD_NEGATIVE_AUDIT_CSV = RESULTS_DIR / "score_robustness" / "audited_hard_negative_audit_summary.csv"
AUDITED_HARD_NEGATIVE_CATEGORIES = {
    "cosmetics",
    "food_contact_substances",
    "human_drugs",
    "animal_drugs",
    "flavoring_agents",
    "fragrances",
    "lipids",
}
FULL_DECISION_CSV = RESULTS_DIR / "pubchem_full_category_decision.csv"

PAPER_EVIDENCE_CSV = PAPER_DIR / "evidence_table.csv"

CHOI_FRAGRANCE_JSON = DATA_MODELS_DIR / "choi_fragrance.json"
CHOI_SURFACTANT_JSON = DATA_MODELS_DIR / "choi_surfactant.json"
CHOI_COSMETIC_ORIGINAL_JSON = OUTPUT_DIR / "experimental_choi_cosmetic.json"
CHOI_COSMETIC_HYBRID_JSON = OUTPUT_DIR / "choi_cosmetic_hybrid_samepair.json"


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))

def load_optional_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    return load_csv_rows(path)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def maybe_float(value: Any) -> float | None:
    if value in (None, "", "None"):
        return None
    return float(value)


def fmt_float(value: float | None, digits: int = 4) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def first_row(rows: list[dict[str, str]], **filters: str) -> dict[str, str]:
    for row in rows:
        if all(str(row.get(key, "")) == str(value) for key, value in filters.items()):
            return row
    raise KeyError(f"Missing row for filters={filters}")


def maybe_first_row(rows: list[dict[str, str]], **filters: str) -> dict[str, str] | None:
    for row in rows:
        if all(str(row.get(key, "")) == str(value) for key, value in filters.items()):
            return row
    return None

def mean(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def sampled_run_counts(per_run_rows: list[dict[str, str]], category: str, regime: str, field: str) -> list[int]:
    return [
        int(row[field])
        for row in per_run_rows
        if row["category"] == category and row["regime"] == regime and row.get(field) not in ("", None)
    ]


def optimization_methods_for_runs(per_run_rows: list[dict[str, str]], category: str, regime: str) -> str:
    methods = sorted(
        {
            row["optimization_method"]
            for row in per_run_rows
            if row["category"] == category
            and row["regime"] == regime
            and row.get("optimization_method") not in ("", None)
        }
    )
    return " | ".join(methods)


def fmt_mean_count(counts: list[int]) -> str:
    value = mean([float(count) for count in counts])
    return str(int(round(value))) if value is not None else ""


def source_artifact_for_stress_regime(regime: str) -> str:
    artifact_paths = [str(SCORE_SUMMARY_CSV)]
    if regime == "audited_hard_negative":
        artifact_paths.append(str(AUDITED_HARD_NEGATIVE_AUDIT_CSV))
    return " | ".join(artifact_paths)



def build_stress_row(
    *,
    category: str,
    display_name: str,
    regime: str,
    regime_label: str,
    negative_definition: str,
    summary_rows: list[dict[str, str]],
    per_run_rows: list[dict[str, str]],
    full_rows: list[dict[str, str]],
    audited_audit_rows: list[dict[str, str]],
) -> dict[str, str] | None:
    summary = maybe_first_row(summary_rows, category=category, regime=regime)
    if summary is None:
        return None
    full = first_row(full_rows, category=category)
    thresholds = [
        float(row["threshold"])
        for row in per_run_rows
        if row["category"] == category and row["regime"] == regime and row.get("threshold") not in ("", None)
    ]
    positive_counts = sampled_run_counts(per_run_rows, category, regime, "positive_count")
    negative_counts = sampled_run_counts(per_run_rows, category, regime, "negative_count")
    notes = (
        "Unrepaired broad hard-negative stress metric; do not compare directly to native/original notebook rows."
        if regime == "property_matched"
        else "Audited broad hard-negative metric after deterministic negative-set audit; keep alongside the unrepaired property-matched row."
    )
    return {
        "category": display_name,
        "benchmark_name": f"{display_name} / {regime_label}",
        "regime_block": f"stress_test_{regime}",
        "positive_definition": "Broad PubChem HID 72 whole-category positives.",
        "negative_definition": negative_definition,
        "optimization_method": optimization_methods_for_runs(per_run_rows, category, regime),
        "positive_count": fmt_mean_count(positive_counts),
        "negative_count": fmt_mean_count(negative_counts),
        "auc": fmt_float(maybe_float(summary.get("auc_mean"))),
        "ks": fmt_float(maybe_float(summary.get("ks_mean"))),
        "balanced_accuracy": fmt_float(maybe_float(summary.get("balanced_accuracy_mean"))),
        "threshold": fmt_float(mean(thresholds)),
        "whole_category_auc_context": fmt_float(maybe_float(full.get("auc"))),
        "source_artifact": source_artifact_for_stress_regime(regime),
        "notes": notes,
    }


def build_stress_rows(
    *,
    category: str,
    display_name: str,
    summary_rows: list[dict[str, str]],
    per_run_rows: list[dict[str, str]],
    full_rows: list[dict[str, str]],
    audited_audit_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    property_row = build_stress_row(
        category=category,
        display_name=display_name,
        regime="property_matched",
        regime_label="unrepaired property-matched stress test",
        negative_definition="Unrepaired property-matched hard negatives drawn from related neighboring categories in robustness validation.",
        summary_rows=summary_rows,
        per_run_rows=per_run_rows,
        full_rows=full_rows,
        audited_audit_rows=audited_audit_rows,
    )
    if property_row is not None:
        rows.append(property_row)
    audited_row = build_stress_row(
        category=category,
        display_name=display_name,
        regime="audited_hard_negative",
        regime_label="audited hard-negative stress test",
        negative_definition="Audited broad hard negatives retained after overlap, dual-use, and near-positive review; not a replacement for the unrepaired property-matched row.",
        summary_rows=summary_rows,
        per_run_rows=per_run_rows,
        full_rows=full_rows,
        audited_audit_rows=audited_audit_rows,
    )
    if audited_row is not None:
        rows.append(audited_row)
    return rows


def build_original_row(
    *,
    category: str,
    benchmark_name: str,
    negative_definition: str,
    optimization_method: str,
    auc: float,
    ks: float | None,
    balanced_accuracy: float | None,
    threshold: float | None,
    source_artifact: str,
    notes: str,
) -> dict[str, str]:
    return {
        "category": category,
        "benchmark_name": benchmark_name,
        "regime_block": "original_submission_benchmark",
        "positive_definition": "Student submission's native positive set.",
        "negative_definition": negative_definition,
        "optimization_method": optimization_method,
        "positive_count": "",
        "negative_count": "",
        "auc": fmt_float(auc),
        "ks": fmt_float(ks),
        "balanced_accuracy": fmt_float(balanced_accuracy),
        "threshold": fmt_float(threshold),
        "whole_category_auc_context": "",
        "source_artifact": source_artifact,
        "notes": notes,
    }


def build_hybrid_row(
    *,
    category: str,
    benchmark_name: str,
    config: dict[str, Any],
    source_artifact: Path,
    notes: str,
) -> dict[str, str]:
    metrics = config.get("metrics", {})
    return {
        "category": category,
        "benchmark_name": benchmark_name,
        "regime_block": "hybrid_pipeline_benchmark",
        "positive_definition": "Same category input pair as the Choi automated pipeline comparison row.",
        "negative_definition": f"Negative source file = {Path(str(config.get('negative_source_csv', ''))).name}; Tanimoto filtering threshold = {config.get('tanimoto_threshold', '')}.",
        "optimization_method": str(config.get("optimization_method", "grid_w_auc_only")),
        "positive_count": "",
        "negative_count": str(metrics.get("negative_count", "")),
        "auc": fmt_float(maybe_float(metrics.get("auc"))),
        "ks": fmt_float(maybe_float(metrics.get("ks"))),
        "balanced_accuracy": fmt_float(maybe_float(metrics.get("balanced_accuracy"))),
        "threshold": fmt_float(maybe_float(config.get("threshold"))),
        "whole_category_auc_context": "",
        "source_artifact": str(source_artifact),
        "notes": notes,
    }


def baseline_regime_label(regime_block: str) -> str:
    if regime_block == "original_submission_benchmark":
        return "original_native_submission"
    if regime_block == "hybrid_pipeline_benchmark":
        return "same_pair_hybrid"
    if regime_block == "stress_test_property_matched":
        return "broad_unrepaired_property_matched"
    if regime_block == "stress_test_audited_hard_negative":
        return "broad_audited_hard_negative"
    return regime_block


def stress_baseline_regime_label(regime: str) -> str:
    if regime == "property_matched":
        return "broad_unrepaired_property_matched"
    if regime == "audited_hard_negative":
        return "broad_audited_hard_negative"
    return f"broad_stress_{regime}"


def build_baseline_matrix(
    rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    per_run_rows: list[dict[str, str]],
    full_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    matrix: list[dict[str, str]] = []
    seen_stress_keys: set[tuple[str, str]] = set()
    seen_audited_categories: set[str] = set()
    for row in rows:
        evidence_regime = baseline_regime_label(row["regime_block"])
        if evidence_regime.startswith("broad_"):
            seen_stress_keys.add((row["category"], evidence_regime))
        if evidence_regime == "broad_audited_hard_negative":
            seen_audited_categories.add(row["category"])
        matrix.append(
            {
                "category": row["category"],
                "evidence_regime": evidence_regime,
                "benchmark_name": row["benchmark_name"],
                "availability": "available",
                "auc": row["auc"],
                "ks": row["ks"],
                "balanced_accuracy": row["balanced_accuracy"],
                "threshold": row["threshold"],
                "positive_count": row["positive_count"],
                "negative_count": row["negative_count"],
                "source_artifact": row["source_artifact"],
                "comparison_note": row["notes"],
            }
        )

    full_by_category = {row["category"]: row for row in full_rows}
    audited_supported_display_names = {
        full_by_category[slug].get("name", slug) if slug in full_by_category else slug
        for slug in AUDITED_HARD_NEGATIVE_CATEGORIES
    }
    for summary in summary_rows:
        category_slug = summary["category"]
        regime = summary["regime"]
        evidence_regime = stress_baseline_regime_label(regime)
        full = full_by_category.get(category_slug)
        display_name = full.get("name", category_slug) if full is not None else category_slug
        if (display_name, evidence_regime) in seen_stress_keys:
            continue
        thresholds = [
            float(row["threshold"])
            for row in per_run_rows
            if row["category"] == category_slug and row["regime"] == regime and row.get("threshold") not in ("", None)
        ]
        positive_counts = sampled_run_counts(per_run_rows, category_slug, regime, "positive_count")
        negative_counts = sampled_run_counts(per_run_rows, category_slug, regime, "negative_count")
        if evidence_regime == "broad_audited_hard_negative":
            seen_audited_categories.add(display_name)
        matrix.append(
            {
                "category": display_name,
                "evidence_regime": evidence_regime,
                "benchmark_name": f"{display_name} / {regime.replace('_', '-')} broad stress baseline",
                "availability": "available",
                "auc": fmt_float(maybe_float(summary.get("auc_mean"))),
                "ks": fmt_float(maybe_float(summary.get("ks_mean"))),
                "balanced_accuracy": fmt_float(maybe_float(summary.get("balanced_accuracy_mean"))),
                "threshold": fmt_float(mean(thresholds)),
                "positive_count": fmt_mean_count(positive_counts),
                "negative_count": fmt_mean_count(negative_counts),
                "source_artifact": source_artifact_for_stress_regime(regime),
                "comparison_note": "Broad stress baseline from score robustness outputs; keep separate from native/original and same-pair hybrid evidence.",
            }
        )
        seen_stress_keys.add((display_name, evidence_regime))

    categories_with_audited_placeholders = {row["category"] for row in matrix} & audited_supported_display_names
    for category in sorted(categories_with_audited_placeholders - seen_audited_categories):
        matrix.append(
            {
                "category": category,
                "evidence_regime": "broad_audited_hard_negative",
                "benchmark_name": f"{category} / audited hard-negative stress test",
                "availability": "not_yet_available",
                "auc": "",
                "ks": "",
                "balanced_accuracy": "",
                "threshold": "",
                "positive_count": "",
                "negative_count": "",
                "source_artifact": str(AUDITED_HARD_NEGATIVE_AUDIT_CSV),
                "comparison_note": "Placeholder baseline slot; populate from audited_hard_negative rows after the audit framework regenerates score robustness outputs.",
            }
        )
    return matrix

def main() -> None:
    summary_rows = load_csv_rows(SCORE_SUMMARY_CSV)
    per_run_rows = load_csv_rows(SCORE_PER_RUN_CSV)
    full_rows = load_csv_rows(FULL_DECISION_CSV)
    audited_audit_rows = load_optional_csv_rows(AUDITED_HARD_NEGATIVE_AUDIT_CSV)

    choi_fragrance = load_json(CHOI_FRAGRANCE_JSON)
    choi_surfactant = load_json(CHOI_SURFACTANT_JSON)
    choi_cosmetic_original = load_json(CHOI_COSMETIC_ORIGINAL_JSON)
    choi_cosmetic_hybrid = load_json(CHOI_COSMETIC_HYBRID_JSON)

    rows: list[dict[str, str]] = []

    rows.append(
        build_original_row(
            category="Pesticides",
            benchmark_name="Pesticides / Kim Nayeon original submission",
            negative_definition="ZINC-derived structurally distant negatives chosen by dynamic Morgan-fingerprint/Tanimoto cutoff crossing-point logic, then balanced by random sampling.",
            optimization_method="submission_specific_ppv_histogram_rule",
            auc=0.9710,
            ks=None,
            balanced_accuracy=None,
            threshold=None,
            source_artifact="김나연_20250786_pesticide.ipynb",
            notes="Original notebook benchmark; PR-AUC 0.9763, accuracy 0.9169, MCC 0.8337 were also reported in the source summary.",
        )
    )
    rows.append(
        build_original_row(
            category="Pesticides",
            benchmark_name="Pesticides / Lee Seoyun original submission",
            negative_definition="Drug-derived negatives after Morgan-fingerprint/Tanimoto filtering with threshold <= 0.24 against agrochemical positives.",
            optimization_method="joint_bayesian_optimization_in_notebook",
            auc=0.9520,
            ks=0.8710,
            balanced_accuracy=0.9350,
            threshold=0.3430,
            source_artifact="이서윤_20251288_Agrochemical.ipynb",
            notes="This is Lee's native benchmark regime, not the broad PubChem stress test.",
        )
    )
    rows.extend(
        build_stress_rows(
            category="pesticides",
            display_name="Pesticides",
            summary_rows=summary_rows,
            per_run_rows=per_run_rows,
            full_rows=full_rows,
            audited_audit_rows=audited_audit_rows,
        )
    )

    rows.append(
        build_original_row(
            category="Cosmetics",
            benchmark_name="Cosmetics / Choi Yebin original automated pipeline",
            negative_definition="Pesticide-derived negatives after Morgan-fingerprint/Tanimoto filtering with threshold < 0.30.",
            optimization_method="grid_w_auc_only",
            auc=0.7080,
            ks=None,
            balanced_accuracy=None,
            threshold=None,
            source_artifact="최예빈_20251266_pesticide.ipynb",
            notes="Original notebook auto-pipeline selected patterns and then optimized only w for AUC.",
        )
    )
    rows.append(
        build_hybrid_row(
            category="Cosmetics",
            benchmark_name="Cosmetics / Choi pipeline + Lee Bayesian optimization (same input pair)",
            config=choi_cosmetic_hybrid,
            source_artifact=CHOI_COSMETIC_HYBRID_JSON,
            notes="Same cosmetics-vs-pesticides input pair as the original Choi notebook, but re-fit with the hybrid multi-objective tuner plus fallback.",
        )
    )
    rows.extend(
        build_stress_rows(
            category="cosmetics",
            display_name="Cosmetics",
            summary_rows=summary_rows,
            per_run_rows=per_run_rows,
            full_rows=full_rows,
            audited_audit_rows=audited_audit_rows,
        )
    )

    rows.append(
        build_original_row(
            category="Fragrances",
            benchmark_name="Fragrances / Choi Yebin original automated pipeline",
            negative_definition="Pesticide-derived negatives after Morgan-fingerprint/Tanimoto filtering with threshold < 0.30.",
            optimization_method="grid_w_auc_only",
            auc=0.8950,
            ks=None,
            balanced_accuracy=None,
            threshold=None,
            source_artifact="최예빈_20251266_pesticide.ipynb",
            notes="Original notebook auto-pipeline with AUC-only optimization of w.",
        )
    )
    rows.append(
        build_hybrid_row(
            category="Fragrances",
            benchmark_name="Fragrances / Choi pipeline + Lee Bayesian optimization (same input pair)",
            config=choi_fragrance,
            source_artifact=CHOI_FRAGRANCE_JSON,
            notes="Same fragrance-vs-pesticides input pair as the original Choi notebook, re-fit with the hybrid multi-objective tuner.",
        )
    )
    rows.extend(
        build_stress_rows(
            category="fragrances",
            display_name="Fragrances",
            summary_rows=summary_rows,
            per_run_rows=per_run_rows,
            full_rows=full_rows,
            audited_audit_rows=audited_audit_rows,
        )
    )

    rows.append(
        build_original_row(
            category="Surfactants",
            benchmark_name="Surfactants / Choi Yebin original automated pipeline",
            negative_definition="Pesticide-derived negatives after Morgan-fingerprint/Tanimoto filtering with threshold < 0.30.",
            optimization_method="grid_w_auc_only",
            auc=0.9660,
            ks=None,
            balanced_accuracy=None,
            threshold=None,
            source_artifact="최예빈_20251266_pesticide.ipynb",
            notes="Original notebook auto-pipeline with AUC-only optimization of w.",
        )
    )
    rows.append(
        build_hybrid_row(
            category="Surfactants",
            benchmark_name="Surfactants / Choi pipeline + Lee Bayesian optimization (same input pair)",
            config=choi_surfactant,
            source_artifact=CHOI_SURFACTANT_JSON,
            notes="Same surfactant-vs-pesticides input pair as the original Choi notebook, re-fit with the hybrid multi-objective tuner.",
        )
    )
    rows.extend(
        build_stress_rows(
            category="surfactants",
            display_name="Surfactants",
            summary_rows=summary_rows,
            per_run_rows=per_run_rows,
            full_rows=full_rows,
            audited_audit_rows=audited_audit_rows,
        )
    )

    fieldnames = list(rows[0].keys())
    for path in (OUT_CSV, OUT_PAPER_CSV):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    OUT_JSON.write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")
    baseline_matrix = build_baseline_matrix(rows, summary_rows, per_run_rows, full_rows)
    baseline_fieldnames = list(baseline_matrix[0].keys())
    with OUT_BASELINE_MATRIX_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=baseline_fieldnames)
        writer.writeheader()
        writer.writerows(baseline_matrix)
    OUT_BASELINE_MATRIX_JSON.write_text(json.dumps(baseline_matrix, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {len(rows)} comparison rows and {len(baseline_matrix)} baseline rows")


if __name__ == "__main__":
    main()
