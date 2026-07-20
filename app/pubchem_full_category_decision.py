from __future__ import annotations

import argparse
import csv
import json
import random
from pathlib import Path
from typing import Any

from build_scoring_models import CHOI_DEFAULT_BAYES_TRIALS, build_choi_model
from pubchem_category_pipeline import (
    CHEMICAL_CLASSES_PATH,
    OUTPUT_DIR,
    RESULTS_DIR,
    fetch_cids_for_hnid,
    fetch_hid72_tree,
    fetch_smiles_for_cids,
    flatten_hid72_tree,
    sample_cids,
    write_smiles_csv,
)

FULL_DECISION_JSON = RESULTS_DIR / "pubchem_full_category_decision.json"
FULL_DECISION_CSV = RESULTS_DIR / "pubchem_full_category_decision.csv"
FULL_DECISION_SUMMARY = RESULTS_DIR / "pubchem_full_category_decision_summary.json"
DECISION_OUTPUT_DIR = OUTPUT_DIR / "full_category_decision"
EXISTING_RELEASE_SLUGS = {"pesticides", "fragrances", "surfactants"}


def _slug(text: str) -> str:
    cleaned = []
    for ch in (text or "").strip().lower():
        cleaned.append(ch if ch.isalnum() else "_")
    slug = "".join(cleaned)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")


def _select_full_targets(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    targets = []
    for row in rows:
        if not row.get("is_leaf"):
            continue
        if not str(row.get("path", "")).startswith(CHEMICAL_CLASSES_PATH + " > "):
            continue
        slug = _slug(str(row["name"]))
        targets.append(
            {
                "name": row["name"],
                "slug": slug,
                "hnid": int(row["hnid"]),
                "compound_count": int(row["compound_count"]),
                "path": row["path"],
                "already_in_existing_release_line": slug in EXISTING_RELEASE_SLUGS,
            }
        )
    targets.sort(key=lambda item: (item["path"], item["name"]))
    return targets


def _decision_label(auc: float | None, balanced_accuracy: float | None, positive_count: int, already_release_line: bool) -> str:
    if positive_count < 25:
        return "insufficient_positive_examples"
    if auc is None or balanced_accuracy is None:
        return "build_failed"
    if already_release_line:
        return "existing_release_line"
    if auc >= 0.90 and balanced_accuracy >= 0.80 and positive_count >= 100:
        return "strong_release_candidate"
    if auc >= 0.80 and balanced_accuracy >= 0.75 and positive_count >= 100:
        return "release_candidate"
    if auc >= 0.70:
        return "experimental_revisit"
    return "reject_for_now"


def run_full_decision(
    max_positive: int,
    max_negative: int,
    min_positive: int,
    bayes_trials: int,
    seed: int,
) -> list[dict[str, Any]]:
    data = fetch_hid72_tree()
    rows = flatten_hid72_tree(data)
    targets = _select_full_targets(rows)

    DECISION_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    inputs_dir = DECISION_OUTPUT_DIR / "inputs"
    models_dir = DECISION_OUTPUT_DIR / "models"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    models_dir.mkdir(parents=True, exist_ok=True)

    smiles_by_slug: dict[str, list[str]] = {}
    for index, target in enumerate(targets):
        cids = fetch_cids_for_hnid(target["hnid"])
        sampled = sample_cids(cids, max_positive, seed + index)
        smiles_by_slug[target["slug"]] = fetch_smiles_for_cids(sampled)

    results: list[dict[str, Any]] = []
    for index, target in enumerate(targets):
        slug = target["slug"]
        positive_smiles = smiles_by_slug.get(slug, [])
        if len(positive_smiles) < min_positive:
            results.append(
                {
                    "model_id": f"full_{slug}",
                    "name": target["name"],
                    "category": slug,
                    "hnid": target["hnid"],
                    "compound_count": target["compound_count"],
                    "positive_count": len(positive_smiles),
                    "negative_count": 0,
                    "auc": None,
                    "ks": None,
                    "balanced_accuracy": None,
                    "objective": None,
                    "threshold": None,
                    "optimization_method": "",
                    "optimization_trials": 0,
                    "decision": _decision_label(None, None, len(positive_smiles), target["already_in_existing_release_line"]),
                    "already_in_existing_release_line": target["already_in_existing_release_line"],
                    "model_path": "",
                }
            )
            continue

        negative_pool: list[str] = []
        for other in targets:
            if other["slug"] == slug:
                continue
            negative_pool.extend(smiles_by_slug.get(other["slug"], []))
        negative_pool = list(dict.fromkeys(negative_pool))
        rng = random.Random(seed + 1000 + index)
        if len(negative_pool) > max_negative:
            negative_pool = rng.sample(negative_pool, max_negative)

        positive_csv = inputs_dir / f"{slug}__positive.csv"
        negative_csv = inputs_dir / f"{slug}__negative.csv"
        write_smiles_csv(positive_csv, positive_smiles)
        write_smiles_csv(negative_csv, negative_pool)

        model_path = models_dir / f"full_{slug}.json"
        build_choi_model(
            positive_csv=positive_csv,
            negative_source_csv=negative_csv,
            model_id=f"full_{slug}",
            label=f"{target['name']} / Full Decision",
            category=slug,
            output_path=model_path,
            bayes_trials=bayes_trials,
            seed=seed + index,
        )
        config = json.loads(model_path.read_text(encoding="utf-8"))
        metrics = config.get("metrics", {})
        auc = float(metrics.get("auc")) if metrics.get("auc") is not None else None
        ks = float(metrics.get("ks")) if metrics.get("ks") is not None else None
        balanced_accuracy = float(metrics.get("balanced_accuracy")) if metrics.get("balanced_accuracy") is not None else None
        objective = float(metrics.get("objective")) if metrics.get("objective") is not None else None
        threshold = float(config.get("threshold")) if config.get("threshold") is not None else None

        results.append(
            {
                "model_id": f"full_{slug}",
                "name": target["name"],
                "category": slug,
                "hnid": target["hnid"],
                "compound_count": target["compound_count"],
                "positive_count": len(positive_smiles),
                "negative_count": len(negative_pool),
                "auc": auc,
                "ks": ks,
                "balanced_accuracy": balanced_accuracy,
                "objective": objective,
                "threshold": threshold,
                "optimization_method": str(config.get("optimization_method", "")),
                "optimization_trials": int(config.get("optimization_trials", 0)),
                "decision": _decision_label(auc, balanced_accuracy, len(positive_smiles), target["already_in_existing_release_line"]),
                "already_in_existing_release_line": target["already_in_existing_release_line"],
                "model_path": str(model_path),
            }
        )
    return results


def save_results(results: list[dict[str, Any]]) -> None:
    FULL_DECISION_JSON.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    with FULL_DECISION_CSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "model_id",
                "name",
                "category",
                "hnid",
                "compound_count",
                "positive_count",
                "negative_count",
                "auc",
                "ks",
                "balanced_accuracy",
                "objective",
                "threshold",
                "optimization_method",
                "optimization_trials",
                "decision",
                "already_in_existing_release_line",
                "model_path",
            ],
        )
        writer.writeheader()
        for row in results:
            payload = dict(row)
            for key in ["auc", "ks", "balanced_accuracy", "objective", "threshold"]:
                if isinstance(payload.get(key), float):
                    payload[key] = f"{payload[key]:.4f}"
            writer.writerow(payload)

    summary = {
        "total_categories": len(results),
        "strong_release_candidate": [row["category"] for row in results if row["decision"] == "strong_release_candidate"],
        "release_candidate": [row["category"] for row in results if row["decision"] == "release_candidate"],
        "experimental_revisit": [row["category"] for row in results if row["decision"] == "experimental_revisit"],
        "reject_for_now": [row["category"] for row in results if row["decision"] == "reject_for_now"],
        "insufficient_positive_examples": [row["category"] for row in results if row["decision"] == "insufficient_positive_examples"],
        "existing_release_line": [row["category"] for row in results if row["decision"] == "existing_release_line"],
    }
    FULL_DECISION_SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build and score the full PubChem HID 72 chemical-class category set for decision support.")
    parser.add_argument("--max-positive", type=int, default=1200)
    parser.add_argument("--max-negative", type=int, default=3000)
    parser.add_argument("--min-positive", type=int, default=25)
    parser.add_argument("--bayes-trials", type=int, default=min(8, CHOI_DEFAULT_BAYES_TRIALS))
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = run_full_decision(
        max_positive=args.max_positive,
        max_negative=args.max_negative,
        min_positive=args.min_positive,
        bayes_trials=args.bayes_trials,
        seed=args.seed,
    )
    save_results(results)
    print(FULL_DECISION_CSV)
    print(FULL_DECISION_SUMMARY)


if __name__ == "__main__":
    main()
