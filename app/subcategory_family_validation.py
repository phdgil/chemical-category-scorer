from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any

from build_evidence_panels import MAX_FAMILIES, PANEL_SPECS, _cluster_entries, _prepare_entries, _read_smiles
from build_scoring_models import build_choi_model
from validate_subtyping_reason import build_entries, build_property_matched_negatives, load_category_smiles, sample_list, write_smiles_csv

APP_DIR = Path(__file__).resolve().parent
EVIDENCE_DIR = APP_DIR / "data" / "evidence_panels"
SOURCE_MAP_JSON = APP_DIR.parent / "results" / "pubchem_named_subtype_sources.json"

def _humanize_family_name(name: str) -> str:
    text = str(name or "prototype family").replace("_", " ").strip()
    return text[:1].upper() + text[1:] if text else "Prototype family"


def _load_source_map() -> dict[str, list[dict[str, Any]]]:
    if not SOURCE_MAP_JSON.exists():
        return {}
    rows = json.loads(SOURCE_MAP_JSON.read_text(encoding="utf-8"))
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row.get("panel_id", "")), []).append(row)
    return grouped


def _panel_family_metadata(panel_id: str) -> dict[str, dict[str, Any]]:
    panel_path = EVIDENCE_DIR / f"{panel_id}.json"
    if not panel_path.exists():
        return {}
    panel = json.loads(panel_path.read_text(encoding="utf-8"))
    source_rows = _load_source_map().get(panel_id, [])
    source_names = [row.get("pubchem_source_name", "") for row in source_rows if row.get("pubchem_source_name")]
    source_hnids = [str(row.get("hnid", "")) for row in source_rows if row.get("hnid")]
    metadata: dict[str, dict[str, Any]] = {}
    for family in panel.get("families", []):
        family_id = str(family.get("family_id", ""))
        family_name = str(family.get("family_name", "prototype_family"))
        metadata[family_id] = {
            "family_name": family_name,
            "display_name": f"{panel_id.replace('_', ' ').title()} prototype: {_humanize_family_name(family_name)}",
            "source_candidates": source_names,
            "source_hnids": source_hnids,
        }
    return metadata

APP_DIR = Path(__file__).resolve().parent
RESULTS_DIR = APP_DIR.parent / "results" / "subcategory_validation"
RUN_DIR = APP_DIR / "output" / "subcategory_validation"
TARGET_PANELS = ["cosmetics", "food_contact_substances", "human_drugs", "animal_drugs"]


def build_panel_families(panel_id: str) -> list[dict[str, Any]]:
    spec = PANEL_SPECS[panel_id]
    smiles_list = _read_smiles(spec["positive_csv"])
    entries = _prepare_entries(smiles_list)
    clusters = _cluster_entries(entries)
    clusters.sort(key=len, reverse=True)
    family_metadata = _panel_family_metadata(panel_id)
    families = []
    for family_index, cluster in enumerate(clusters[:MAX_FAMILIES], start=1):
        smiles = [entries[idx]["smiles"] for idx in cluster]
        family_id = f"family_{family_index:02d}"
        meta = family_metadata.get(family_id, {})
        family_name = str(meta.get("family_name", f"prototype_family_{family_index}"))
        display_name = str(meta.get("display_name", f"{panel_id.replace('_', ' ').title()} prototype: {_humanize_family_name(family_name)}"))
        families.append(
            {
                "panel_id": panel_id,
                "family_id": family_id,
                "family_name": family_name,
                "display_name": display_name,
                "source_candidates": meta.get("source_candidates", []),
                "source_hnids": meta.get("source_hnids", []),
                "member_count": len(smiles),
                "smiles": list(dict.fromkeys(smiles)),
            }
        )
    return families


def family_negative_regimes(
    panel_id: str,
    family_smiles: list[str],
    categories: dict[str, list[str]],
    category_entries: dict[str, list[dict[str, Any]]],
    max_negative: int,
    seed: int,
) -> dict[str, list[str]]:
    parent_smiles = list(dict.fromkeys(categories[panel_id]))
    family_set = set(family_smiles)
    sibling_smiles = [smiles for smiles in parent_smiles if smiles not in family_set]

    external_smiles = []
    external_entries = []
    for slug, smiles_list in categories.items():
        if slug == panel_id:
            continue
        external_smiles.extend(smiles_list)
        external_entries.extend(category_entries.get(slug, []))
    external_smiles = list(dict.fromkeys(external_smiles))

    all_other_smiles = list(dict.fromkeys(sibling_smiles + external_smiles))
    positive_entries = build_entries(family_smiles)
    property_matched = build_property_matched_negatives(positive_entries, external_entries + category_entries.get(panel_id, []), max_negative)

    return {
        "all_other_random": sample_list(all_other_smiles, max_negative, seed + 11),
        "sibling_hard": sample_list(sibling_smiles, max_negative, seed + 22),
        "external_other": sample_list(external_smiles, max_negative, seed + 33),
        "property_matched": property_matched,
    }


def build_family_model(
    panel_id: str,
    family_id: str,
    regime: str,
    seed: int,
    positive_smiles: list[str],
    negative_smiles: list[str],
    bayes_trials: int,
) -> dict[str, Any]:
    run_dir = RUN_DIR / panel_id / family_id / regime / f"seed_{seed}"
    positive_csv = run_dir / "positive.csv"
    negative_csv = run_dir / "negative.csv"
    model_path = run_dir / "model.json"
    write_smiles_csv(positive_csv, positive_smiles)
    write_smiles_csv(negative_csv, negative_smiles)
    build_choi_model(
        positive_csv=positive_csv,
        negative_source_csv=negative_csv,
        model_id=f"subcat_{panel_id}_{family_id}_{regime}_{seed}",
        label=f"{panel_id} {family_id} {regime} seed {seed}",
        category=f"{panel_id}_{family_id}",
        output_path=model_path,
        bayes_trials=bayes_trials,
        seed=seed,
    )
    return json.loads(model_path.read_text(encoding="utf-8"))


def summarize(values: list[float]) -> dict[str, float]:
    return {
        "mean": round(statistics.mean(values), 4),
        "sd": round(statistics.stdev(values), 4) if len(values) > 1 else 0.0,
        "min": round(min(values), 4),
        "max": round(max(values), 4),
    }


def run_validation(seeds: list[int], max_negative: int, bayes_trials: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    categories = load_category_smiles()
    category_entries = {slug: build_entries(smiles) for slug, smiles in categories.items()}
    per_run_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    for panel_id in TARGET_PANELS:
        families = build_panel_families(panel_id)
        for family in families:
            regime_rows: dict[str, list[dict[str, Any]]] = {regime: [] for regime in ["all_other_random", "sibling_hard", "external_other", "property_matched"]}
            for seed in seeds:
                regimes = family_negative_regimes(
                    panel_id=panel_id,
                    family_smiles=family["smiles"],
                    categories=categories,
                    category_entries=category_entries,
                    max_negative=max_negative,
                    seed=seed,
                )
                for regime, negative_smiles in regimes.items():
                    if not negative_smiles:
                        continue
                    config = build_family_model(
                        panel_id=panel_id,
                        family_id=family["family_id"],
                        regime=regime,
                        seed=seed,
                        positive_smiles=family["smiles"],
                        negative_smiles=negative_smiles,
                        bayes_trials=bayes_trials,
                    )
                    metrics = config.get("metrics", {})
                    row = {
                        "panel_id": panel_id,
                        "family_id": family["family_id"],
                        "family_name": family["family_name"],
                        "display_name": family["display_name"],
                        "member_count": family["member_count"],
                        "seed": seed,
                        "regime": regime,
                        "negative_count": int(metrics.get("negative_count", 0)),
                        "auc": round(float(metrics.get("auc", 0.0)), 4),
                        "ks": round(float(metrics.get("ks", 0.0)), 4),
                        "balanced_accuracy": round(float(metrics.get("balanced_accuracy", 0.0)), 4),
                        "overlap": round(float(metrics.get("overlap", 0.0)), 4),
                        "threshold": round(float(config.get("threshold", 0.5)), 4),
                        "optimization_method": str(config.get("optimization_method", "")),
                    }

                    per_run_rows.append(row)
                    regime_rows[regime].append(row)

            for regime, rows in regime_rows.items():
                if not rows:
                    continue
                auc_stats = summarize([float(row["auc"]) for row in rows])
                ba_stats = summarize([float(row["balanced_accuracy"]) for row in rows])
                summary_rows.append(
                    {
                        "panel_id": panel_id,
                        "family_id": family["family_id"],
                        "family_name": family["family_name"],
                        "display_name": family["display_name"],
                        "member_count": family["member_count"],
                        "regime": regime,
                        "runs": len(rows),
                        "auc_mean": auc_stats["mean"],
                        "auc_sd": auc_stats["sd"],
                        "auc_min": auc_stats["min"],
                        "auc_max": auc_stats["max"],
                        "balanced_accuracy_mean": ba_stats["mean"],
                        "balanced_accuracy_sd": ba_stats["sd"],
                    }
                )
    return per_run_rows, summary_rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate performance of evidence-panel subcategories (prototype families).")
    parser.add_argument("--seeds", type=int, nargs="+", default=[11, 23, 37])
    parser.add_argument("--max-negative", type=int, default=3000)
    parser.add_argument("--bayes-trials", type=int, default=3)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    per_run_rows, summary_rows = run_validation(seeds=args.seeds, max_negative=args.max_negative, bayes_trials=args.bayes_trials)
    write_csv(
        RESULTS_DIR / "subcategory_per_run.csv",
        per_run_rows,
        [
            "panel_id",
            "family_id",
            "family_name",
            "display_name",
            "member_count",
            "seed",
            "regime",
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
        RESULTS_DIR / "subcategory_summary.csv",
        summary_rows,
        [
            "panel_id",
            "family_id",
            "family_name",
            "display_name",
            "member_count",
            "regime",
            "runs",
            "auc_mean",
            "auc_sd",
            "auc_min",
            "auc_max",
            "balanced_accuracy_mean",
            "balanced_accuracy_sd",
        ],
    )
    (RESULTS_DIR / "subcategory_summary.json").write_text(
        json.dumps({"seeds": args.seeds, "summary": summary_rows}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(RESULTS_DIR)


if __name__ == "__main__":
    main()
