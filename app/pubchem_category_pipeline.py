from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from algorithm_score_engine import MODELS_DIR, list_models, refresh_model_registry
from build_scoring_models import CHOI_DEFAULT_BAYES_TRIALS, build_choi_model


APP_DIR = Path(__file__).resolve().parent
DATA_DIR = APP_DIR / "data"
OUTPUT_DIR = APP_DIR / "output" / "pubchem_pipeline"
RESULTS_DIR = APP_DIR.parent / "results"
CATALOG_JSON_PATH = RESULTS_DIR / "pubchem_hid72_catalog.json"
CATALOG_CSV_PATH = RESULTS_DIR / "pubchem_hid72_catalog.csv"
TARGETS_JSON_PATH = RESULTS_DIR / "pubchem_build_targets.json"
TARGETS_CSV_PATH = RESULTS_DIR / "pubchem_build_targets.csv"
SCREEN_JSON_PATH = RESULTS_DIR / "pubchem_category_screen.json"
SCREEN_CSV_PATH = RESULTS_DIR / "pubchem_category_screen.csv"
REGISTRY_PATH = RESULTS_DIR / "model_registry.json"
RELEASE_CONFIG_PATH = DATA_DIR / "app_release_config.json"
SUMMARY_PATH = RESULTS_DIR / "algorithm_project_summary.json"
RAW_TREE_PATH = RESULTS_DIR / "pubchem_hid72_raw.json"

TREE_URL = "https://pubchem.ncbi.nlm.nih.gov/classification/cgi/classifications.fcgi?format=json&hid=72"
CID_URL_TEMPLATE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/classification/hnid/{hnid}/cids/TXT"
PROPERTY_URL_TEMPLATE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/ConnectivitySMILES/CSV"


ROOT_NAME = "PubChem Compound TOC"
CHEMICAL_CLASSES_PATH = f"{ROOT_NAME} > Chemical and Physical Properties > Chemical Classes"

DEFAULT_MIN_COUNT = 200
DEFAULT_MAX_POSITIVE = 3500
DEFAULT_MAX_NEGATIVE = 5000
DEFAULT_RELEASE_AUC = 0.80
DEFAULT_SEED = 20260710
CID_BATCH_SIZE = 200
PROPERTY_BATCH_SIZE = 100
REQUEST_SLEEP_SECONDS = 0.25
DEFAULT_BAYES_TRIALS = min(8, CHOI_DEFAULT_BAYES_TRIALS)

EXISTING_RELEASE_CATEGORY_SLUGS = {"pesticides", "fragrances", "surfactants"}



def _slug(text: str) -> str:
    cleaned = []
    for ch in (text or "").strip().lower():
        if ch.isalnum():
            cleaned.append(ch)
        else:
            cleaned.append("_")
    slug = "".join(cleaned)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")


def _request_json(url: str, retries: int = 4) -> Any:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as response:
                return json.load(response)
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep((attempt + 1) * 1.0)
    raise RuntimeError(f"Failed to fetch JSON: {url}")


def _request_text(url: str, retries: int = 4) -> str:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as response:
                return response.read().decode("utf-8")
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep((attempt + 1) * 1.0)
    raise RuntimeError(f"Failed to fetch text: {url}")


def fetch_hid72_tree() -> dict[str, Any]:
    data = _request_json(TREE_URL)
    RAW_TREE_PATH.parent.mkdir(parents=True, exist_ok=True)
    RAW_TREE_PATH.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    return data


def flatten_hid72_tree(data: dict[str, Any]) -> list[dict[str, Any]]:
    hierarchy = data["Hierarchies"]["Hierarchy"][0]
    nodes = hierarchy["Node"]
    by_id = {node["NodeID"]: node for node in nodes}

    def node_name(node: dict[str, Any]) -> str:
        name = node.get("Information", {}).get("Name", {})
        if isinstance(name, dict):
            markup = name.get("StringWithMarkup", {})
            if isinstance(markup, dict):
                return str(markup.get("String", ""))
            if isinstance(markup, list):
                return " ".join(str(item.get("String", "")) for item in markup if isinstance(item, dict)).strip()
        return str(name or "")

    rows: list[dict[str, Any]] = []

    def walk(node_id: str, depth: int, path_parts: list[str]) -> None:
        node = by_id[node_id]
        info = node.get("Information", {})
        child_ids = [child_id for child_id in info.get("ChildID", []) if child_id in by_id]
        counts = info.get("Counts", [])
        count = int(counts[0].get("Count", 0)) if counts else 0
        name = node_name(node)
        path = path_parts + [name]
        rows.append(
            {
                "node_id": node_id,
                "hnid": int(info.get("HNID")) if info.get("HNID") is not None else None,
                "name": name,
                "depth": depth,
                "compound_count": count,
                "child_count": len(child_ids),
                "is_leaf": len(child_ids) == 0,
                "path": " > ".join(path),
                "parent_path": " > ".join(path_parts),
            }
        )
        for child_id in child_ids:
            walk(child_id, depth + 1, path)

    root_info = hierarchy.get("Information", {})
    for child_id in root_info.get("ChildID", []):
        if child_id in by_id:
            walk(child_id, 1, [ROOT_NAME])
    return rows


def save_catalog(rows: list[dict[str, Any]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    CATALOG_JSON_PATH.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    with CATALOG_CSV_PATH.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["node_id", "hnid", "name", "depth", "compound_count", "child_count", "is_leaf", "parent_path", "path"],
        )
        writer.writeheader()
        writer.writerows(rows)


def select_build_targets(rows: list[dict[str, Any]], min_count: int) -> list[dict[str, Any]]:
    targets = []
    for row in rows:
        if not str(row["path"]).startswith(CHEMICAL_CLASSES_PATH + " > "):
            continue
        if not row["is_leaf"]:
            continue
        if int(row["compound_count"]) < min_count:
            continue
        targets.append(
            {
                "name": row["name"],
                "slug": _slug(row["name"]),
                "hnid": row["hnid"],
                "compound_count": row["compound_count"],
                "path": row["path"],
                "release_eligible": _slug(row["name"]) not in EXISTING_RELEASE_CATEGORY_SLUGS,
            }
        )
    targets.sort(key=lambda item: (-int(item["compound_count"]), item["name"]))
    return targets


def save_targets(targets: list[dict[str, Any]]) -> None:
    TARGETS_JSON_PATH.write_text(json.dumps(targets, ensure_ascii=False, indent=2), encoding="utf-8")
    with TARGETS_CSV_PATH.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "slug", "hnid", "compound_count", "path", "release_eligible"])
        writer.writeheader()
        writer.writerows(targets)


def fetch_cids_for_hnid(hnid: int) -> list[int]:
    text = _request_text(CID_URL_TEMPLATE.format(hnid=hnid))
    cids = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            cids.append(int(line))
        except ValueError:
            continue
    return cids


def fetch_smiles_for_cids(cids: list[int]) -> list[str]:
    smiles: list[str] = []
    for start in range(0, len(cids), PROPERTY_BATCH_SIZE):
        batch = cids[start : start + PROPERTY_BATCH_SIZE]
        cid_text = ",".join(str(cid) for cid in batch)
        url = PROPERTY_URL_TEMPLATE.format(cids=cid_text)
        text = _request_text(url)
        reader = csv.DictReader(text.splitlines())
        for row in reader:
            value = (row.get("ConnectivitySMILES") or row.get("CanonicalSMILES") or row.get("IsomericSMILES") or "").strip()
            if value:
                smiles.append(value)
        time.sleep(REQUEST_SLEEP_SECONDS)
    deduped = list(dict.fromkeys(smiles))
    return deduped


def sample_cids(cids: list[int], limit: int, seed: int) -> list[int]:
    unique = list(dict.fromkeys(cids))
    if len(unique) <= limit:
        return unique
    rng = random.Random(seed)
    return sorted(rng.sample(unique, limit))


def write_smiles_csv(path: Path, smiles: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["SMILES"])
        writer.writeheader()
        for value in smiles:
            writer.writerow({"SMILES": value})


def build_target_models(
    targets: list[dict[str, Any]],
    max_positive: int,
    max_negative: int,
    release_auc: float,
    seed: int,
    bayes_trials: int,
) -> list[dict[str, Any]]:

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    inputs_dir = OUTPUT_DIR / "inputs"
    models_dir = OUTPUT_DIR / "models"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    models_dir.mkdir(parents=True, exist_ok=True)
    for existing_path in MODELS_DIR.glob("pubchem_*.json"):
        existing_path.unlink(missing_ok=True)

    sampled_smiles: dict[str, list[str]] = {}
    sampled_counts: dict[str, int] = {}

    for index, target in enumerate(targets):
        cids = fetch_cids_for_hnid(int(target["hnid"]))
        sampled_cids = sample_cids(cids, max_positive, seed + index)
        smiles = fetch_smiles_for_cids(sampled_cids)
        sampled_smiles[target["slug"]] = smiles
        sampled_counts[target["slug"]] = len(smiles)

    results: list[dict[str, Any]] = []
    for index, target in enumerate(targets):
        slug = target["slug"]
        positive_smiles = sampled_smiles.get(slug, [])
        if len(positive_smiles) < 100:
            results.append(
                {
                    "model_id": f"pubchem_{slug}",
                    "name": target["name"],
                    "category": slug,
                    "hnid": target["hnid"],
                    "positive_count": len(positive_smiles),
                    "negative_count": 0,
                    "auc": None,
                    "ks": None,
                    "balanced_accuracy": None,
                    "objective": None,
                    "threshold": None,
                    "optimization_method": "",
                    "optimization_trials": 0,
                    "status": "skipped_too_small",
                    "promoted": False,
                    "release_eligible": target["release_eligible"],
                    "model_path": "",
                }

            )
            continue

        negative_pool: list[str] = []
        for other in targets:
            if other["slug"] == slug:
                continue
            negative_pool.extend(sampled_smiles.get(other["slug"], []))
        negative_pool = list(dict.fromkeys(negative_pool))
        rng = random.Random(seed + 1000 + index)
        if len(negative_pool) > max_negative:
            negative_pool = rng.sample(negative_pool, max_negative)
        if len(negative_pool) < 100:
            results.append(
                {
                    "model_id": f"pubchem_{slug}",
                    "name": target["name"],
                    "category": slug,
                    "hnid": target["hnid"],
                    "positive_count": len(positive_smiles),
                    "negative_count": len(negative_pool),
                    "auc": None,
                    "ks": None,
                    "balanced_accuracy": None,
                    "objective": None,
                    "threshold": None,
                    "optimization_method": "",
                    "optimization_trials": 0,
                    "status": "skipped_negative_too_small",
                    "promoted": False,
                    "release_eligible": target["release_eligible"],
                    "model_path": "",
                }
            )
            continue

        positive_csv = inputs_dir / f"{slug}__positive.csv"
        negative_csv = inputs_dir / f"{slug}__negative_source.csv"
        write_smiles_csv(positive_csv, positive_smiles)
        write_smiles_csv(negative_csv, negative_pool)

        experimental_model_path = models_dir / f"pubchem_{slug}.json"
        build_choi_model(
            positive_csv=positive_csv,
            negative_source_csv=negative_csv,
            model_id=f"pubchem_{slug}",
            label=f"{target['name']} / PubChem Auto",
            category=slug,
            output_path=experimental_model_path,
            bayes_trials=bayes_trials,
            seed=seed + index,
        )

        config = json.loads(experimental_model_path.read_text(encoding="utf-8"))
        metrics = config.get("metrics", {})
        auc_value = metrics.get("auc")
        auc = float(auc_value) if auc_value is not None else None
        ks = float(metrics.get("ks", 0.0)) if metrics.get("ks") is not None else None
        balanced_accuracy = float(metrics.get("balanced_accuracy", 0.0)) if metrics.get("balanced_accuracy") is not None else None
        objective = float(metrics.get("objective", 0.0)) if metrics.get("objective") is not None else None
        threshold = float(config.get("threshold", 0.0))
        optimization_method = str(config.get("optimization_method", ""))
        optimization_trials = int(config.get("optimization_trials", 0))

        promoted = bool(auc is not None and auc >= release_auc and target["release_eligible"])
        status = "release" if promoted else ("experimental" if auc is not None else "failed")
        if promoted:
            shutil.copy2(experimental_model_path, MODELS_DIR / experimental_model_path.name)

        results.append(
            {
                "model_id": f"pubchem_{slug}",
                "name": target["name"],
                "category": slug,
                "hnid": target["hnid"],
                "positive_count": len(positive_smiles),
                "negative_count": len(negative_pool),
                "auc": auc,
                "ks": ks,
                "balanced_accuracy": balanced_accuracy,
                "objective": objective,
                "threshold": threshold,
                "optimization_method": optimization_method,
                "optimization_trials": optimization_trials,
                "status": status,
                "promoted": promoted,
                "release_eligible": target["release_eligible"],
                "model_path": str(experimental_model_path),
            }
        )
    return results


def save_screen_results(results: list[dict[str, Any]]) -> None:
    SCREEN_JSON_PATH.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    with SCREEN_CSV_PATH.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "model_id",
                "name",
                "category",
                "hnid",
                "positive_count",
                "negative_count",
                "auc",
                "ks",
                "balanced_accuracy",
                "objective",
                "threshold",
                "optimization_method",
                "optimization_trials",
                "status",
                "promoted",
                "release_eligible",
                "model_path",
            ],

        )
        writer.writeheader()
        for row in results:
            payload = dict(row)
            if isinstance(payload.get("auc"), float):
                payload["auc"] = f"{payload['auc']:.4f}"
            if isinstance(payload.get("ks"), float):
                payload["ks"] = f"{payload['ks']:.4f}"
            if isinstance(payload.get("balanced_accuracy"), float):
                payload["balanced_accuracy"] = f"{payload['balanced_accuracy']:.4f}"
            if isinstance(payload.get("objective"), float):
                payload["objective"] = f"{payload['objective']:.4f}"
            if isinstance(payload.get("threshold"), float):
                payload["threshold"] = f"{payload['threshold']:.4f}"
            writer.writerow(payload)



def refresh_registry() -> list[dict[str, Any]]:
    refresh_model_registry()
    models = list_models()
    REGISTRY_PATH.write_text(json.dumps(models, ensure_ascii=False, indent=2), encoding="utf-8")
    return models


def update_release_config(models: list[dict[str, Any]]) -> None:
    if RELEASE_CONFIG_PATH.exists():
        config = json.loads(RELEASE_CONFIG_PATH.read_text(encoding="utf-8"))
    else:
        config = {"app_name": "Chemical Category Scorer"}
    model_ids = [model["model_id"] for model in models]
    config["status"] = "v0.5 pubchem-bayesian validated-model release"
    config["default_model"] = "kim_pesticide"
    config["supports_all_model_ranking"] = True
    config["supports_pubchem_class_builder"] = True
    config["available_models"] = model_ids
    config["notes"] = [
        "Built-in scorers now combine validated student baselines with validated PubChem auto-built category models.",
        "The PubChem HID 72 category catalog is saved under results/pubchem_hid72_catalog.csv and .json.",
        "Lee Seoyun-style Bayesian optimization is applied when auto-building PubChem category scorers.",
        "Experimental PubChem-category outputs remain under app/output/pubchem_pipeline/models and results/pubchem_category_screen.csv."
    ]
    RELEASE_CONFIG_PATH.write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")



def update_project_summary(models: list[dict[str, Any]], build_results: list[dict[str, Any]]) -> None:
    if SUMMARY_PATH.exists():
        summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    else:
        summary = {"project": "algorithm_paper_app"}

    integrated_models = []
    for model in models:
        integrated_models.append(
            {
                "model_id": model["model_id"],
                "student": model.get("source_student", "PubChem Auto"),
                "category": model["category"],
                "method": model.get("description", "PubChem auto-built scorer"),
            }
        )
    promoted = [row for row in build_results if row.get("promoted")]
    experimental = [row for row in build_results if row.get("status") == "experimental"]
    summary["integrated_models"] = integrated_models
    summary["pubchem_execution"] = {
        "catalog_path": str(CATALOG_CSV_PATH),
        "target_path": str(TARGETS_CSV_PATH),
        "screen_path": str(SCREEN_CSV_PATH),
        "promoted_model_ids": [row["model_id"] for row in promoted],
        "experimental_model_ids": [row["model_id"] for row in experimental],
    }
    summary["new_capabilities"] = [
        "Desktop app can switch between validated built-in scoring models.",
        "Desktop app can score all built-in models and export comparative batch CSVs.",
        "PubChem HID 72 categories are cataloged automatically before execution.",
        "PubChem chemical-class targets can be auto-built, screened, and promoted into the app when AUC passes the release gate.",
        "Lee Seoyun-style Bayesian optimization now tunes PubChem auto-built scorers before release promotion.",
        "Desktop app can build new Choi-style scorers directly from PubChem classification CSV exports."
    ]
    SUMMARY_PATH.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


def run_discovery(min_count: int) -> list[dict[str, Any]]:
    data = fetch_hid72_tree()
    rows = flatten_hid72_tree(data)
    save_catalog(rows)
    targets = select_build_targets(rows, min_count=min_count)
    save_targets(targets)
    return targets


def run_execution(
    min_count: int,
    max_positive: int,
    max_negative: int,
    release_auc: float,
    seed: int,
    bayes_trials: int,
) -> list[dict[str, Any]]:
    targets = run_discovery(min_count=min_count)
    results = build_target_models(
        targets=targets,
        max_positive=max_positive,
        max_negative=max_negative,
        release_auc=release_auc,
        seed=seed,
        bayes_trials=bayes_trials,
    )

    save_screen_results(results)
    models = refresh_registry()
    update_release_config(models)
    update_project_summary(models, results)
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Catalog PubChem HID 72 categories and build auto-scored category models.")
    parser.add_argument("command", choices=["discover", "execute", "all"])
    parser.add_argument("--min-count", type=int, default=DEFAULT_MIN_COUNT)
    parser.add_argument("--max-positive", type=int, default=DEFAULT_MAX_POSITIVE)
    parser.add_argument("--max-negative", type=int, default=DEFAULT_MAX_NEGATIVE)
    parser.add_argument("--release-auc", type=float, default=DEFAULT_RELEASE_AUC)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--bayes-trials", type=int, default=DEFAULT_BAYES_TRIALS)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "discover":
        targets = run_discovery(min_count=args.min_count)
        print(f"catalog_rows={sum(1 for _ in json.loads(CATALOG_JSON_PATH.read_text(encoding='utf-8')))}")
        print(f"target_count={len(targets)}")
        return
    if args.command in {"execute", "all"}:
        results = run_execution(
            min_count=args.min_count,
            max_positive=args.max_positive,
            max_negative=args.max_negative,
            release_auc=args.release_auc,
            seed=args.seed,
            bayes_trials=args.bayes_trials,
        )
        promoted = [row["model_id"] for row in results if row.get("promoted")]
        print(f"target_count={len(results)}")
        print(f"promoted={','.join(promoted)}")
        return
    raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
