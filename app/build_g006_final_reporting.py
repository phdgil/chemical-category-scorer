from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
RESULTS_DIR = ROOT_DIR / "results"
PAPER_DIR = ROOT_DIR / "paper"
MODELS_DIR = APP_DIR / "data" / "models"
APP_DATA_DIR = APP_DIR / "data"

FULL_DECISION_CSV = RESULTS_DIR / "pubchem_full_category_decision.csv"
ROBUSTNESS_SUMMARY_CSV = RESULTS_DIR / "score_robustness" / "score_robustness_summary.csv"
G003_SUMMARY_CSV = RESULTS_DIR / "g003_broad_category_improvement" / "g003_broad_category_improvement_summary.csv"
G004_SUMMARY_CSV = RESULTS_DIR / "g004_broad_drug_improvement" / "g004_broad_category_improvement_summary.csv"
G005_SUMMARY_CSV = RESULTS_DIR / "g005_remaining_category_improvement" / "g005_remaining_category_improvement_summary.csv"
OFFICIAL_DRUG_SUBTYPE_SUMMARY_CSV = RESULTS_DIR / "official_drug_subtype_validation" / "official_drug_subtype_partial_summary.csv"
SUBCATEGORY_SUMMARY_CSV = RESULTS_DIR / "subcategory_validation" / "subcategory_summary.csv"

FINAL_CSV = RESULTS_DIR / "final_scoring_performance_table.csv"
FINAL_JSON = RESULTS_DIR / "final_scoring_performance_table.json"
PAPER_FINAL_CSV = PAPER_DIR / "final_scoring_performance_table.csv"
PROJECT_SUMMARY_JSON = RESULTS_DIR / "algorithm_project_summary.json"
MODEL_REGISTRY_JSON = RESULTS_DIR / "model_registry.json"
APP_RELEASE_CONFIG_JSON = APP_DATA_DIR / "app_release_config.json"

CATEGORY_ORDER = [
    "pfas",
    "surfactants",
    "solvents",
    "endocrine_disruptors",
    "pesticides",
    "fragrances",
    "food_additives",
    "flavoring_agents",
    "lipids",
    "cosmetics",
    "food_contact_substances",
    "human_drugs",
    "animal_drugs",
    "polymers",
    "uvcb",
]

APP_MODEL_IDS = [
    "choi_fragrance",
    "choi_surfactant",
    "kim_pesticide",
    "lee_pesticide",
    "pubchem_endocrine_disruptors",
    "pubchem_flavoring_agents",
    "pubchem_food_additives",
    "pubchem_lipids",
    "pubchem_pfas",
    "pubchem_solvents",
]

MODEL_CATEGORY_ALIASES = {
    "fragrance": "fragrances",
    "surfactant": "surfactants",
    "pesticide": "pesticides",
}

DISPLAY_OVERRIDES = {
    "pfas": "PFAS",
    "uvcb": "UVCB",
}

G003_CATEGORIES = {"cosmetics", "food_contact_substances"}
G004_CATEGORIES = {"human_drugs", "animal_drugs"}
G005_CATEGORIES = {
    "pesticides",
    "fragrances",
    "food_additives",
    "flavoring_agents",
    "lipids",
    "solvents",
    "endocrine_disruptors",
    "surfactants",
}

SCORER_REFERENCES = [
    {
        "output_id": "pesticides__kim_nayeon_original_scorer",
        "output_name": "Pesticides / Kim Nayeon original scorer",
        "parent_category": "Pesticides",
        "reported_auc": "0.9710",
        "reported_balanced_accuracy": "",
        "source": "kim_original_submission",
        "caveat": "Original notebook benchmark (ROC-AUC 0.9710, PR-AUC 0.9763, Accuracy 0.9169, MCC 0.8337); not comparable to current broad PubChem stress tests.",
    },
    {
        "output_id": "pesticides__lee_seoyun_original_scorer",
        "output_name": "Pesticides / Lee Seoyun original scorer",
        "parent_category": "Pesticides",
        "reported_auc": "0.9520",
        "reported_balanced_accuracy": "0.9350",
        "source": "lee_original_submission",
        "caveat": "Original Lee agrochemical/drug-negative benchmark with Tanimoto <= 0.24, KS 0.871, threshold 0.343; not comparable to current broad PubChem stress tests.",
    },
    {
        "output_id": "cosmetics__choi_yebin_original_scorer",
        "output_name": "Cosmetics / Choi Yebin original scorer",
        "parent_category": "Cosmetics",
        "reported_auc": "0.7080",
        "reported_balanced_accuracy": "",
        "source": "choi_original_submission",
        "caveat": "Original Choi automated scorer used pesticide-derived negatives; retained only as a native-scorer reference beside the weak broad cosmetics target.",
    },
    {
        "output_id": "fragrances__choi_yebin_original_scorer",
        "output_name": "Fragrances / Choi Yebin original scorer",
        "parent_category": "Fragrances",
        "reported_auc": "0.8950",
        "reported_balanced_accuracy": "",
        "source": "choi_original_submission",
        "caveat": "Original Choi automated scorer used pesticide-derived negatives; retained only as a native-scorer reference beside current broad-category caveats.",
    },
    {
        "output_id": "surfactants__choi_yebin_original_scorer",
        "output_name": "Surfactants / Choi Yebin original scorer",
        "parent_category": "Surfactants",
        "reported_auc": "0.9660",
        "reported_balanced_accuracy": "",
        "source": "choi_original_submission",
        "caveat": "Original Choi automated scorer used pesticide-derived negatives; retained only as a native-scorer reference beside current broad-category caveats.",
    },
]

SUPPORTING_PANEL_CAVEAT = "Broad parent status remains caveated by latest bounded broad-category evidence."

SUPPORTING_PANEL_MIN_AUC = 0.8
SUPPORTING_PANEL_MIN_BALANCED_ACCURACY = 0.78



def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: Any) -> str:
    if value in (None, ""):
        return ""
    return f"{float(value):.4f}".rstrip("0").rstrip(".")


def slugify(text: str) -> str:
    chars: list[str] = []
    for ch in (text or "").strip().lower():
        chars.append(ch if ch.isalnum() else "_")
    slug = "".join(chars)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")


def display_name(slug: str, full_by_category: dict[str, dict[str, str]]) -> str:
    if slug in DISPLAY_OVERRIDES:
        return DISPLAY_OVERRIDES[slug]
    row = full_by_category.get(slug)
    if row:
        return row["name"]
    return slug.replace("_", " ").title()


def by_category_and_regime(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    return {(row["category"], row.get("regime", "")): row for row in rows}


def best_variant(rows: list[dict[str, str]], category: str) -> dict[str, str] | None:
    candidates = [row for row in rows if row["category"] == category]
    if not candidates:
        return None
    return max(candidates, key=lambda row: (float(row.get("auc_mean") or 0), float(row.get("balanced_accuracy_mean") or 0)))


def baseline_variant(rows: list[dict[str, str]], category: str) -> dict[str, str] | None:
    for row in rows:
        if row["category"] == category and row["variant"].startswith("baseline"):
            return row
    return None


def bounded_source(category: str) -> str:
    if category in G003_CATEGORIES:
        return "results/g003_broad_category_improvement/g003_broad_category_improvement_summary.csv"
    if category in G004_CATEGORIES:
        return "results/g004_broad_drug_improvement/g004_broad_category_improvement_summary.csv"
    if category in G005_CATEGORIES:
        return "results/g005_remaining_category_improvement/g005_remaining_category_improvement_summary.csv"
    return "results/score_robustness/score_robustness_summary.csv"


def select_current_bounded_row(category: str, summaries: dict[str, list[dict[str, str]]], robustness: dict[tuple[str, str], dict[str, str]]) -> tuple[dict[str, str] | None, dict[str, str] | None, str]:
    if category in G003_CATEGORIES:
        return best_variant(summaries["g003"], category), baseline_variant(summaries["g003"], category), "audited_hard_negative"
    if category in G004_CATEGORIES:
        return best_variant(summaries["g004"], category), baseline_variant(summaries["g004"], category), "audited_hard_negative"
    if category in G005_CATEGORIES:
        best = best_variant(summaries["g005"], category)
        base = baseline_variant(summaries["g005"], category)
        return best, base, (best or base or {}).get("eval_regime") or (best or base or {}).get("regime") or "property_matched"
    return robustness.get((category, "property_matched")) or robustness.get((category, "audited_hard_negative")), None, "property_matched"


def policy_for(category: str, auc: float | None, ba: float | None, source: str, caveat_basis: str) -> str:
    if category in {"polymers", "uvcb"}:
        return "not_prioritized_broad_category"
    if category in G003_CATEGORIES | G004_CATEGORIES:
        return "preserved_weak_broad_category_with_supporting_panels"
    if source == "results/score_robustness/score_robustness_summary.csv" and auc is not None and auc >= 0.9 and (ba or 0) >= 0.85:
        return "robust_built_in_score"
    if auc is not None and auc >= 0.75 and (ba or 0) >= 0.7:
        return "bounded_candidate_do_not_promote_without_multi_seed_confirmation"
    if "positive_delta" in caveat_basis:
        return "experimental_broad_category_attempt_positive_delta"
    return "preserved_weak_broad_category"


def caveat_for(category: str, current: dict[str, str] | None, baseline: dict[str, str] | None, source: str) -> tuple[str, str]:
    if current is None:
        return "No current bounded stress row found; preserve only as an inventory target.", "missing_current_bounded_row"
    runs = current.get("runs", "1")
    variant = current.get("variant", current.get("regime", "property_matched"))
    delta_auc = float(current.get("delta_auc_vs_baseline_mean") or 0)
    delta_ba = float(current.get("delta_ba_vs_baseline_mean") or 0)
    status = "positive_delta" if delta_auc > 0 or delta_ba > 0 else "no_confirmed_improvement"
    if source == "results/score_robustness/score_robustness_summary.csv":
        return f"Current score-robustness row is one bounded run for {current.get('regime', 'property_matched')}; use as release evidence only with the single-seed caveat.", "robustness_single_seed"
    return (
        f"Latest bounded attempt uses variant={variant}, runs={runs}, delta_auc={delta_auc:.4f}, delta_balanced_accuracy={delta_ba:.4f}; this is a bounded/single-seed experimental result and must not be treated as final proof.",
        status,
    )

def supporting_panel_rows(latest_by_display: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    row_template = {key: "" for key in latest_by_display[next(iter(latest_by_display))]}

    for source_row in read_csv(OFFICIAL_DRUG_SUBTYPE_SUMMARY_CSV):
        if source_row.get("regime") != "property_matched":
            continue
        if float(source_row.get("auc_mean") or 0) < SUPPORTING_PANEL_MIN_AUC:
            continue
        if float(source_row.get("balanced_accuracy_mean") or 0) < SUPPORTING_PANEL_MIN_BALANCED_ACCURACY:
            continue
        parent = latest_by_display.get(display_name(source_row["panel_id"], {}))
        if parent is None:
            continue
        subtype_name = source_row["subtype_name"]
        row = dict(row_template)
        row.update(
            {
                "output_id": f"{source_row['panel_id']}__{source_row['subtype_hnid']}__{slugify(subtype_name)}",
                "output_name": f"{parent['output_name']} / {subtype_name}",
                "output_level": "subcategory",
                "parent_category": parent["output_name"],
                "reported_positive_count": "",
                "reported_auc": fmt(source_row["auc_mean"]),
                "reported_balanced_accuracy": fmt(source_row["balanced_accuracy_mean"]),
                "validation_regime": source_row["regime"],
                "recommended_final_policy": "official_subtype_panel",
                "source": "results/official_drug_subtype_validation/official_drug_subtype_partial_summary.csv",
                "parent_whole_category_positive_count": parent["parent_whole_category_positive_count"],
                "parent_whole_category_auc": parent["parent_whole_category_auc"],
                "parent_whole_category_balanced_accuracy": parent["parent_whole_category_balanced_accuracy"],
                "parent_property_matched_auc": parent["parent_property_matched_auc"],
                "parent_property_matched_balanced_accuracy": parent["parent_property_matched_balanced_accuracy"],
                "latest_bounded_broad_status": parent["latest_bounded_broad_status"],
                "latest_bounded_broad_caveat": parent["latest_bounded_broad_caveat"],
                "latest_bounded_source": parent["latest_bounded_source"],
                "sample_size_risk": "subtype_count_not_reported_in_summary",
                "key_interpretation": f"Whole-category score is weak, but official FDA subtype labels recover useful signal for this named subtype. {SUPPORTING_PANEL_CAVEAT}",
            }
        )
        rows.append(row)

    for source_row in read_csv(SUBCATEGORY_SUMMARY_CSV):
        if source_row.get("regime") != "property_matched":
            continue
        if float(source_row.get("auc_mean") or 0) < SUPPORTING_PANEL_MIN_AUC:
            continue
        if float(source_row.get("balanced_accuracy_mean") or 0) < SUPPORTING_PANEL_MIN_BALANCED_ACCURACY:
            continue
        parent = latest_by_display.get(display_name(source_row["panel_id"], {}))
        if parent is None:
            continue
        family_name = source_row["family_name"]
        row = dict(row_template)
        row.update(
            {
                "output_id": f"{source_row['panel_id']}__{source_row['family_id']}__{slugify(family_name)}",
                "output_name": f"{parent['output_name']} / {family_name.replace('_', ' ').title()}",
                "output_level": "subcategory",
                "parent_category": parent["output_name"],
                "reported_positive_count": source_row["member_count"],
                "reported_auc": fmt(source_row["auc_mean"]),
                "reported_balanced_accuracy": fmt(source_row["balanced_accuracy_mean"]),
                "validation_regime": source_row["regime"],
                "recommended_final_policy": "evidence_panel",
                "source": "results/subcategory_validation/subcategory_summary.csv",
                "parent_whole_category_positive_count": parent["parent_whole_category_positive_count"],
                "parent_whole_category_auc": parent["parent_whole_category_auc"],
                "parent_whole_category_balanced_accuracy": parent["parent_whole_category_balanced_accuracy"],
                "parent_property_matched_auc": parent["parent_property_matched_auc"],
                "parent_property_matched_balanced_accuracy": parent["parent_property_matched_balanced_accuracy"],
                "latest_bounded_broad_status": parent["latest_bounded_broad_status"],
                "latest_bounded_broad_caveat": parent["latest_bounded_broad_caveat"],
                "latest_bounded_source": parent["latest_bounded_source"],
                "sample_size_risk": "adequate_size" if int(source_row["member_count"]) >= 10 else "small_sample_or_sparse_category",
                "key_interpretation": f"Whole-category score is weak, but subtype-family scores are strong; report as evidence panel with prototype-family support. {SUPPORTING_PANEL_CAVEAT}",
            }
        )
        rows.append(row)

    return sorted(
        rows,
        key=lambda row: (
            CATEGORY_ORDER.index(slugify(row["parent_category"])) if slugify(row["parent_category"]) in CATEGORY_ORDER else len(CATEGORY_ORDER),
            row["recommended_final_policy"],
            -float(row["reported_auc"] or 0),
            row["output_id"],
        ),
    )


def build_final_rows() -> list[dict[str, str]]:
    full_rows = read_csv(FULL_DECISION_CSV)
    full_by_category = {row["category"]: row for row in full_rows}
    robustness = by_category_and_regime(read_csv(ROBUSTNESS_SUMMARY_CSV))
    summaries = {
        "g003": read_csv(G003_SUMMARY_CSV),
        "g004": read_csv(G004_SUMMARY_CSV),
        "g005": read_csv(G005_SUMMARY_CSV),
    }

    rows: list[dict[str, str]] = []
    latest_by_display: dict[str, dict[str, str]] = {}
    for category in CATEGORY_ORDER:
        full = full_by_category[category]
        current, baseline, regime = select_current_bounded_row(category, summaries, robustness)
        source = bounded_source(category)
        caveat, caveat_basis = caveat_for(category, current, baseline, source)
        auc = float(current.get("auc_mean") or current.get("auc") or full.get("auc") or 0) if current else float(full.get("auc") or 0)
        ba = float(current.get("balanced_accuracy_mean") or current.get("balanced_accuracy") or full.get("balanced_accuracy") or 0) if current else float(full.get("balanced_accuracy") or 0)
        policy = policy_for(category, auc, ba, source, caveat_basis)
        name = display_name(category, full_by_category)
        whole_auc = fmt(full["auc"])
        whole_ba = fmt(full["balanced_accuracy"])
        current_auc = fmt(auc)
        current_ba = fmt(ba)
        prop = robustness.get((category, "property_matched"))
        row = {
            "output_id": category,
            "output_name": name,
            "output_level": "category",
            "parent_category": "",
            "reported_positive_count": full["positive_count"],
            "reported_auc": current_auc,
            "reported_balanced_accuracy": current_ba,
            "validation_regime": regime,
            "recommended_final_policy": policy,
            "source": "broad_category_current_evidence",
            "parent_whole_category_positive_count": full["positive_count"],
            "parent_whole_category_auc": whole_auc,
            "parent_whole_category_balanced_accuracy": whole_ba,
            "parent_property_matched_auc": fmt(prop.get("auc_mean")) if prop else "",
            "parent_property_matched_balanced_accuracy": fmt(prop.get("balanced_accuracy_mean")) if prop else "",
            "latest_bounded_broad_status": caveat_basis,
            "latest_bounded_broad_caveat": caveat,
            "latest_bounded_source": source,
            "sample_size_risk": "adequate_size" if int(full["positive_count"]) >= 100 else "small_sample_or_sparse_category",
            "key_interpretation": "Broad category remains the reporting target; scorer status is limited by the latest bounded stress evidence and explicit caveats.",
        }
        rows.append(row)
        latest_by_display[name] = row

    for ref in SCORER_REFERENCES:
        parent = latest_by_display[ref["parent_category"]]
        scorer = dict(parent)
        scorer.update(
            {
                "output_id": ref["output_id"],
                "output_name": ref["output_name"],
                "output_level": "scoring_function",
                "parent_category": ref["parent_category"],
                "reported_positive_count": "",
                "reported_auc": ref["reported_auc"],
                "reported_balanced_accuracy": ref["reported_balanced_accuracy"],
                "validation_regime": "original_submission_negative_set",
                "recommended_final_policy": "submission_scorer_reference",
                "source": ref["source"],
                "latest_bounded_broad_status": parent["latest_bounded_broad_status"],
                "latest_bounded_broad_caveat": parent["latest_bounded_broad_caveat"],
                "latest_bounded_source": parent["latest_bounded_source"],
                "sample_size_risk": "submission_dataset_not_recast_here",
                "key_interpretation": ref["caveat"],
            }
        )
        rows.append(scorer)

    rows.extend(supporting_panel_rows(latest_by_display))

    return rows


def build_model_registry(final_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    by_output_id = {row["output_id"]: row for row in final_rows if row["output_level"] == "category"}
    registry: list[dict[str, Any]] = []
    for path in sorted(MODELS_DIR.glob("*.json")):
        model = json.loads(path.read_text(encoding="utf-8"))
        model_id = model["model_id"]
        if model_id not in APP_MODEL_IDS:
            continue
        category = model.get("category", "")
        reporting_category = MODEL_CATEGORY_ALIASES.get(category, category)
        status_row = by_output_id.get(reporting_category)
        registry.append(
            {
                "model_id": model_id,
                "label": model.get("label", model_id),
                "category": category,
                "reporting_category": reporting_category,
                "threshold": model.get("threshold"),
                "model_type": model.get("model_type"),
                "source_student": model.get("source_student", ""),
                "description": model.get("description", ""),
                "built_in_app_model": True,
                "release_policy": status_row.get("recommended_final_policy", "student_baseline_reference") if status_row else "student_baseline_reference",
                "broad_category_caveat": status_row.get("latest_bounded_broad_caveat", "Student/native model; broad PubChem status is tracked separately.") if status_row else "Student/native model; broad PubChem status is tracked separately.",
            }
        )
    return sorted(registry, key=lambda row: APP_MODEL_IDS.index(row["model_id"]))


def build_project_summary(registry: list[dict[str, Any]], final_rows: list[dict[str, str]]) -> dict[str, Any]:
    categories = [row for row in final_rows if row["output_level"] == "category"]
    return {
        "project": "algorithm_paper_app",
        "created_for": "D:/DSWU/2026_기말고사",
        "focus": "algorithm paper plus desktop app",
        "g006_refresh": {
            "status": "refreshed_from_post_g005_evidence",
            "caveat": "G003-G005 improvement outputs are bounded one-seed experimental attempts unless explicitly noted; they update reporting status but do not silently promote candidate models.",
            "source_artifacts": [
                "results/score_robustness/score_robustness_summary.csv",
                "results/g003_broad_category_improvement/g003_broad_category_improvement_summary.csv",
                "results/g004_broad_drug_improvement/g004_broad_category_improvement_summary.csv",
                "results/g005_remaining_category_improvement/g005_remaining_category_improvement_summary.csv",
                "results/official_drug_subtype_validation/",
                "results/subcategory_validation/",
            ],
        },
        "selected_baseline": {
            "student": "김나연",
            "student_id": "20250786",
            "model_id": "kim_pesticide",
            "source_notebook": "D:/DSWU/2026_기말고사/컴퓨터알고리즘/김나연/github_repo/김나연_20250786_pesticide.ipynb",
            "metrics": {"roc_auc": 0.971, "pr_auc": 0.9763, "accuracy": 0.9169, "mcc": 0.8337, "decision_threshold": 0.434},
        },
        "integrated_models": [
            {"model_id": row["model_id"], "student": row.get("source_student", ""), "category": row.get("category", ""), "method": row.get("description", ""), "release_policy": row.get("release_policy", "")} for row in registry
        ],
        "broad_category_reporting": [
            {
                "category": row["output_id"],
                "name": row["output_name"],
                "policy": row["recommended_final_policy"],
                "latest_bounded_status": row["latest_bounded_broad_status"],
                "latest_bounded_caveat": row["latest_bounded_broad_caveat"],
                "source": row["latest_bounded_source"],
            }
            for row in categories
        ],
        "app_model_policy": {
            "built_in_model_ids": APP_MODEL_IDS,
            "experimental_candidate_scope": "G003-G005 candidate model artifacts remain experimental tracking outputs and are not promoted into app/data/models by this reporting refresh.",
        },
    }


def build_app_release_config(registry: list[dict[str, Any]], final_rows: list[dict[str, str]]) -> dict[str, Any]:
    categories = [row for row in final_rows if row["output_level"] == "category"]
    weak = [row["output_id"] for row in categories if "weak" in row["recommended_final_policy"]]
    experimental = [row["output_id"] for row in categories if "experimental" in row["recommended_final_policy"]]
    robust = [row["output_id"] for row in categories if "robust" in row["recommended_final_policy"]]
    bounded_candidate = [row["output_id"] for row in categories if "bounded_candidate" in row["recommended_final_policy"]]
    not_prioritized = [row["output_id"] for row in categories if "not_prioritized" in row["recommended_final_policy"]]
    return {
        "app_name": "Chemical Category Scorer",
        "status": "v0.6 post-G005 reporting refresh; built-in models unchanged, experimental broad-category tracking explicit",
        "default_model": "kim_pesticide",
        "supports_all_model_ranking": True,
        "supports_pubchem_class_builder": True,
        "available_models": [row["model_id"] for row in registry],
        "built_in_model_policy": "Regenerated from current app/data/models only; G003-G005 experimental candidates are tracked but not promoted into the built-in app set.",
        "post_g005_broad_category_state": {
            "robust_or_release_grade_with_caveats": robust,
            "preserved_weak_categories": weak,
            "experimental_positive_delta_categories": experimental,
            "bounded_candidate_categories": bounded_candidate,
            "not_prioritized_categories": not_prioritized,
            "single_seed_caveat": "G003-G005 rows are bounded one-seed attempts; report as current evidence, not final proof.",
        },
        "tracking_artifacts": [
            "results/final_scoring_performance_table.csv",
            "results/final_scoring_performance_table.json",
            "results/score_robustness/",
            "results/g003_broad_category_improvement/",
            "results/g004_broad_drug_improvement/",
            "results/g005_remaining_category_improvement/",
            "results/official_drug_subtype_validation/",
            "results/subcategory_validation/",
        ],
        "notes": [
            "Broad categories remain the main reporting targets; subtype and original student scorer rows are supporting references only.",
            "Kim/Lee pesticide and Choi native scorers remain visible where they explain original student performance regimes.",
            "Weak broad categories are preserved in reporting and app status instead of being dropped or silently replaced by subtype panels.",
        ],
    }


def main() -> None:
    final_rows = build_final_rows()
    write_csv(FINAL_CSV, final_rows)
    write_csv(PAPER_FINAL_CSV, final_rows)
    FINAL_JSON.write_text(json.dumps(final_rows, indent=2, ensure_ascii=False), encoding="utf-8")

    registry = build_model_registry(final_rows)
    MODEL_REGISTRY_JSON.write_text(json.dumps(registry, indent=2, ensure_ascii=False), encoding="utf-8")

    PROJECT_SUMMARY_JSON.write_text(json.dumps(build_project_summary(registry, final_rows), indent=2, ensure_ascii=False), encoding="utf-8")
    APP_RELEASE_CONFIG_JSON.write_text(json.dumps(build_app_release_config(registry, final_rows), indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {len(final_rows)} final reporting rows and {len(registry)} registry entries")


if __name__ == "__main__":
    main()
