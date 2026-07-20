from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.request
from pathlib import Path
from typing import Any

from build_scoring_models import CHOI_DEFAULT_BAYES_TRIALS, build_choi_model
from pubchem_category_pipeline import PROPERTY_BATCH_SIZE, PROPERTY_URL_TEMPLATE, REQUEST_SLEEP_SECONDS, fetch_cids_for_hnid
from validate_subtyping_reason import build_entries, build_property_matched_negatives, load_category_smiles, sample_list, write_smiles_csv

APP_DIR = Path(__file__).resolve().parent
RESULTS_DIR = APP_DIR.parent / "results" / "official_drug_subtype_validation"
RUN_DIR = APP_DIR / "output" / "official_drug_subtype_validation"
CATALOG_CSV = RESULTS_DIR / "official_drug_subtype_catalog.csv"
CATALOG_JSON = RESULTS_DIR / "official_drug_subtype_catalog.json"
PER_RUN_CSV = RESULTS_DIR / "official_drug_subtype_per_run.csv"
SUMMARY_CSV = RESULTS_DIR / "official_drug_subtype_summary.csv"
SUMMARY_JSON = RESULTS_DIR / "official_drug_subtype_summary.json"

CLASSIFICATION2_URL = "https://pubchem.ncbi.nlm.nih.gov/classification_2/classification_2.fcgi?hid={hid}&start=root&format=json&depth={depth}"

DRUG_SOURCE_SPECS: dict[str, dict[str, Any]] = {
    "human_drugs": {
        "source_hid": 116,
        "source_name": "FDA Drug Type and Pharmacologic Classification",
        "hid72_parent_hnid": 12950172,
        "fetch_depth": 3,
        "candidate_depths": {2, 3},
        "leaf_only": True,
        "external_categories": [
            "animal_drugs",
            "endocrine_disruptors",
            "food_additives",
            "flavoring_agents",
            "solvents",
            "pesticides",
        ],
    },
    "animal_drugs": {
        "source_hid": 136,
        "source_name": "ATCvet Classification",
        "hid72_parent_hnid": 12950173,
        "fetch_depth": 4,
        "candidate_depths": {1, 2, 3},
        "leaf_only": False,
        "external_categories": [
            "human_drugs",
            "endocrine_disruptors",
            "food_additives",
            "flavoring_agents",
            "solvents",
            "pesticides",
        ],
    },
}


def _request_json(url: str, retries: int = 4) -> Any:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as response:
                text = response.read().decode("utf-8", "replace")
                return json.loads(text)
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep((attempt + 1) * 1.0)
    raise RuntimeError(f"Failed to fetch JSON: {url}")


def _request_text(url: str, retries: int = 4) -> str:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as response:
                return response.read().decode("utf-8", "replace")
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep((attempt + 1) * 1.0)
    raise RuntimeError(f"Failed to fetch text: {url}")


def _node_name(info: dict[str, Any]) -> str:
    name = info.get("Name", {})
    if isinstance(name, dict):
        markup = name.get("StringWithMarkup", {})
        if isinstance(markup, dict):
            return str(markup.get("String", ""))
        if isinstance(markup, list):
            return " ".join(str(item.get("String", "")) for item in markup if isinstance(item, dict)).strip()
    return str(name or "")


def _slug(text: str) -> str:
    chars: list[str] = []
    for ch in (text or "").strip().lower():
        chars.append(ch if ch.isalnum() else "_")
    slug = "".join(chars)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")


def _dedupe(values: list[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def fetch_source_tree(hid: int, depth: int) -> dict[str, Any]:
    return _request_json(CLASSIFICATION2_URL.format(hid=hid, depth=depth))["Hierarchies"]["Hierarchy"][0]


def flatten_source_tree(tree: dict[str, Any]) -> list[dict[str, Any]]:
    nodes = tree.get("Node", [])
    by_id = {node["NodeID"]: node for node in nodes}
    children: dict[str, list[str]] = {}
    for node in nodes:
        for parent_id in node.get("ParentID", []):
            children.setdefault(parent_id, []).append(node["NodeID"])

    def depth_and_path(node_id: str) -> tuple[int, str]:
        parts: list[str] = []
        current = by_id[node_id]
        while True:
            parts.append(_node_name(current.get("Information", {})))
            parents = current.get("ParentID", [])
            if not parents or parents[0] == "root":
                break
            current = by_id[parents[0]]
        parts.reverse()
        return len(parts), " > ".join(parts)

    rows: list[dict[str, Any]] = []
    for node in nodes:
        info = node.get("Information", {})
        counts = info.get("Counts", [])
        count = int(counts[0].get("Count", 0)) if counts else 0
        depth, path = depth_and_path(node["NodeID"])
        rows.append(
            {
                "node_id": node["NodeID"],
                "hnid": int(info.get("HNID")) if info.get("HNID") is not None else None,
                "name": _node_name(info),
                "depth": depth,
                "path": path,
                "source_compound_count": count,
                "child_count": len(children.get(node["NodeID"], [])),
                "is_leaf": len(children.get(node["NodeID"], [])) == 0,
            }
        )
    rows.sort(key=lambda row: (row["depth"], row["path"]))
    return rows


def fetch_cid_smiles_map(cids: list[int]) -> dict[int, str]:
    cid_to_smiles: dict[int, str] = {}
    for start in range(0, len(cids), PROPERTY_BATCH_SIZE):
        batch = cids[start : start + PROPERTY_BATCH_SIZE]
        cid_text = ",".join(str(cid) for cid in batch)
        text = _request_text(PROPERTY_URL_TEMPLATE.format(cids=cid_text))
        reader = csv.DictReader(text.splitlines())
        for row in reader:
            cid_value = str(row.get("CID", "")).strip()
            smiles = (row.get("ConnectivitySMILES") or row.get("CanonicalSMILES") or row.get("IsomericSMILES") or "").strip()
            if cid_value and smiles:
                try:
                    cid_to_smiles[int(cid_value)] = smiles
                except ValueError:
                    continue
        time.sleep(REQUEST_SLEEP_SECONDS)
    return cid_to_smiles


def summarize(values: list[float]) -> dict[str, float]:
    if not values:
        return {"mean": 0.0, "sd": 0.0, "min": 0.0, "max": 0.0}
    if len(values) == 1:
        value = round(values[0], 4)
        return {"mean": value, "sd": 0.0, "min": value, "max": value}
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    return {
        "mean": round(mean, 4),
        "sd": round(variance ** 0.5, 4),
        "min": round(min(values), 4),
        "max": round(max(values), 4),
    }


def build_catalog(min_overlap_smiles: int) -> tuple[list[dict[str, Any]], dict[str, list[str]], dict[int, list[str]]]:
    rows: list[dict[str, Any]] = []
    parent_smiles_by_panel: dict[str, list[str]] = {}
    subtype_smiles_cache: dict[int, list[str]] = {}

    for panel_id, spec in DRUG_SOURCE_SPECS.items():
        parent_cids = list(dict.fromkeys(fetch_cids_for_hnid(int(spec["hid72_parent_hnid"]))))
        parent_cid_set = set(parent_cids)
        cid_to_smiles = fetch_cid_smiles_map(parent_cids)
        parent_smiles = _dedupe([cid_to_smiles[cid] for cid in parent_cids if cid in cid_to_smiles])
        parent_smiles_by_panel[panel_id] = parent_smiles

        tree = fetch_source_tree(int(spec["source_hid"]), int(spec["fetch_depth"]))
        nodes = flatten_source_tree(tree)
        for node in nodes:
            if node["hnid"] is None or int(node["source_compound_count"]) <= 0:
                continue
            try:
                subtype_cids = list(dict.fromkeys(fetch_cids_for_hnid(int(node["hnid"]))))
            except Exception:
                continue
            overlap_cids = [cid for cid in subtype_cids if cid in parent_cid_set]
            overlap_smiles = _dedupe([cid_to_smiles[cid] for cid in overlap_cids if cid in cid_to_smiles])
            subtype_smiles_cache[int(node["hnid"])] = overlap_smiles

            rows.append(
                {
                    "panel_id": panel_id,
                    "source_name": spec["source_name"],
                    "source_hid": int(spec["source_hid"]),
                    "hid72_parent_hnid": int(spec["hid72_parent_hnid"]),
                    "subtype_hnid": int(node["hnid"]),
                    "subtype_name": str(node["name"]),
                    "subtype_slug": _slug(str(node["name"])),
                    "subtype_depth": int(node["depth"]),
                    "subtype_path": str(node["path"]),
                    "is_leaf": bool(node["is_leaf"]),
                    "source_compound_count": int(node["source_compound_count"]),
                    "hid72_parent_cid_count": len(parent_cids),
                    "hid72_overlap_cid_count": len(overlap_cids),
                    "hid72_overlap_unique_smiles_count": len(overlap_smiles),
                    "representative_smiles": overlap_smiles[0] if overlap_smiles else "",
                    "selected_for_validation": (
                        int(node["depth"]) in spec["candidate_depths"]
                        and (not spec["leaf_only"] or bool(node["is_leaf"]))
                        and len(overlap_smiles) >= min_overlap_smiles
                    ),
                }
            )

    rows.sort(
        key=lambda row: (
            row["panel_id"],
            not bool(row["selected_for_validation"]),
            -int(row["hid72_overlap_unique_smiles_count"]),
            row["subtype_path"],
        )
    )
    return rows, parent_smiles_by_panel, subtype_smiles_cache


def validation_categories(parent_smiles_by_panel: dict[str, list[str]]) -> dict[str, list[str]]:
    categories = load_category_smiles()
    for panel_id, smiles in parent_smiles_by_panel.items():
        categories[panel_id] = smiles
    return categories


def negative_regimes(
    panel_id: str,
    subtype_smiles: list[str],
    parent_smiles: list[str],
    categories: dict[str, list[str]],
    max_negative: int,
    seed: int,
) -> dict[str, list[str]]:
    positive_set = set(subtype_smiles)
    sibling_smiles = [smiles for smiles in parent_smiles if smiles not in positive_set]
    external_smiles: list[str] = []
    for slug in DRUG_SOURCE_SPECS[panel_id]["external_categories"]:
        external_smiles.extend(categories.get(slug, []))
    external_smiles = _dedupe(external_smiles)
    all_other_smiles = _dedupe(sibling_smiles + external_smiles)
    positive_entries = build_entries(subtype_smiles)
    candidate_entries = build_entries(_dedupe(sibling_smiles + external_smiles))
    property_matched = build_property_matched_negatives(positive_entries, candidate_entries, max_negative)
    return {
        "all_other_random": sample_list(all_other_smiles, max_negative, seed + 11),
        "sibling_hard": sample_list(sibling_smiles, max_negative, seed + 22),
        "external_other": sample_list(external_smiles, max_negative, seed + 33),
        "property_matched": property_matched,
    }


def build_subtype_model(
    panel_id: str,
    subtype_hnid: int,
    subtype_slug: str,
    regime: str,
    seed: int,
    positive_smiles: list[str],
    negative_smiles: list[str],
    bayes_trials: int,
) -> dict[str, Any]:
    run_dir = RUN_DIR / panel_id / f"{subtype_hnid}_{subtype_slug}" / regime / f"seed_{seed}"
    positive_csv = run_dir / "positive.csv"
    negative_csv = run_dir / "negative.csv"
    model_path = run_dir / "model.json"
    write_smiles_csv(positive_csv, positive_smiles)
    write_smiles_csv(negative_csv, negative_smiles)
    build_choi_model(
        positive_csv=positive_csv,
        negative_source_csv=negative_csv,
        model_id=f"official_{panel_id}_{subtype_hnid}_{regime}_{seed}",
        label=f"{panel_id} {subtype_hnid} {regime} seed {seed}",
        category=f"{panel_id}_{subtype_hnid}",
        output_path=model_path,
        bayes_trials=bayes_trials,
        seed=seed,
    )
    return json.loads(model_path.read_text(encoding="utf-8"))


def run_validation(
    catalog_rows: list[dict[str, Any]],
    parent_smiles_by_panel: dict[str, list[str]],
    subtype_smiles_cache: dict[int, list[str]],
    categories: dict[str, list[str]],
    seeds: list[int],
    max_negative: int,
    bayes_trials: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    per_run_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    for row in [item for item in catalog_rows if item["selected_for_validation"]]:
        panel_id = str(row["panel_id"])
        subtype_hnid = int(row["subtype_hnid"])
        subtype_smiles = subtype_smiles_cache[subtype_hnid]
        parent_smiles = parent_smiles_by_panel[panel_id]
        regime_rows: dict[str, list[dict[str, Any]]] = {key: [] for key in ["all_other_random", "sibling_hard", "external_other", "property_matched"]}

        for seed in seeds:
            regimes = negative_regimes(
                panel_id=panel_id,
                subtype_smiles=subtype_smiles,
                parent_smiles=parent_smiles,
                categories=categories,
                max_negative=max_negative,
                seed=seed,
            )
            for regime, negative_smiles in regimes.items():
                if not negative_smiles:
                    continue
                config = build_subtype_model(
                    panel_id=panel_id,
                    subtype_hnid=subtype_hnid,
                    subtype_slug=str(row["subtype_slug"]),
                    regime=regime,
                    seed=seed,
                    positive_smiles=subtype_smiles,
                    negative_smiles=negative_smiles,
                    bayes_trials=bayes_trials,
                )
                metrics = config.get("metrics", {})
                payload = {
                    "panel_id": panel_id,
                    "source_name": row["source_name"],
                    "source_hid": row["source_hid"],
                    "hid72_parent_hnid": row["hid72_parent_hnid"],
                    "subtype_hnid": subtype_hnid,
                    "subtype_name": row["subtype_name"],
                    "subtype_depth": row["subtype_depth"],
                    "seed": seed,
                    "regime": regime,
                    "positive_unique_smiles_count": len(subtype_smiles),
                    "negative_count": int(metrics.get("negative_count", 0)),
                    "auc": round(float(metrics.get("auc", 0.0)), 4),
                    "ks": round(float(metrics.get("ks", 0.0)), 4),
                    "balanced_accuracy": round(float(metrics.get("balanced_accuracy", 0.0)), 4),
                    "overlap": round(float(metrics.get("overlap", 0.0)), 4),
                    "threshold": round(float(config.get("threshold", 0.5)), 4),
                    "optimization_method": str(config.get("optimization_method", "")),
                }
                per_run_rows.append(payload)
                regime_rows[regime].append(payload)

        for regime, rows in regime_rows.items():
            if not rows:
                continue
            auc_stats = summarize([float(item["auc"]) for item in rows])
            ba_stats = summarize([float(item["balanced_accuracy"]) for item in rows])
            summary_rows.append(
                {
                    "panel_id": panel_id,
                    "source_name": row["source_name"],
                    "source_hid": row["source_hid"],
                    "hid72_parent_hnid": row["hid72_parent_hnid"],
                    "subtype_hnid": subtype_hnid,
                    "subtype_name": row["subtype_name"],
                    "subtype_depth": row["subtype_depth"],
                    "positive_unique_smiles_count": len(subtype_smiles),
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

    summary_rows.sort(key=lambda item: (item["panel_id"], item["subtype_name"], item["regime"]))
    return per_run_rows, summary_rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate official human/animal drug PubChem subtype hierarchies against HID72 category membership.")
    parser.add_argument("--seeds", type=int, nargs="+", default=[11])
    parser.add_argument("--max-negative", type=int, default=3000)
    parser.add_argument("--bayes-trials", type=int, default=min(2, CHOI_DEFAULT_BAYES_TRIALS))
    parser.add_argument("--min-overlap-smiles", type=int, default=25)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    catalog_rows, parent_smiles_by_panel, subtype_smiles_cache = build_catalog(min_overlap_smiles=args.min_overlap_smiles)
    categories = validation_categories(parent_smiles_by_panel)

    write_csv(
        CATALOG_CSV,
        catalog_rows,
        [
            "panel_id",
            "source_name",
            "source_hid",
            "hid72_parent_hnid",
            "subtype_hnid",
            "subtype_name",
            "subtype_slug",
            "subtype_depth",
            "subtype_path",
            "is_leaf",
            "source_compound_count",
            "hid72_parent_cid_count",
            "hid72_overlap_cid_count",
            "hid72_overlap_unique_smiles_count",
            "representative_smiles",
            "selected_for_validation",
        ],
    )
    CATALOG_JSON.write_text(json.dumps(catalog_rows, ensure_ascii=False, indent=2), encoding="utf-8")

    per_run_rows, summary_rows = run_validation(
        catalog_rows=catalog_rows,
        parent_smiles_by_panel=parent_smiles_by_panel,
        subtype_smiles_cache=subtype_smiles_cache,
        categories=categories,
        seeds=args.seeds,
        max_negative=args.max_negative,
        bayes_trials=args.bayes_trials,
    )

    write_csv(
        PER_RUN_CSV,
        per_run_rows,
        [
            "panel_id",
            "source_name",
            "source_hid",
            "hid72_parent_hnid",
            "subtype_hnid",
            "subtype_name",
            "subtype_depth",
            "seed",
            "regime",
            "positive_unique_smiles_count",
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
        SUMMARY_CSV,
        summary_rows,
        [
            "panel_id",
            "source_name",
            "source_hid",
            "hid72_parent_hnid",
            "subtype_hnid",
            "subtype_name",
            "subtype_depth",
            "positive_unique_smiles_count",
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
    SUMMARY_JSON.write_text(
        json.dumps(
            {
                "seeds": args.seeds,
                "min_overlap_smiles": args.min_overlap_smiles,
                "catalog": catalog_rows,
                "summary": summary_rows,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(CATALOG_CSV)
    print(SUMMARY_CSV)


if __name__ == "__main__":
    main()
