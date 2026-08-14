from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from sklearn.metrics import roc_auc_score, roc_curve

from app.build_scoring_models import _choi_prepare_rows, _choi_two_scores, build_choi_model
from app.qed_inspired_validation import ROOT, category_positive_smiles

SOURCE_CATEGORIES = (
    "cosmetics",
    "endocrine_disruptors",
    "food_additives",
    "food_contact_substances",
    "pesticides",
    "solvents",
    "surfactants",
)
TARGET_CATEGORIES = ("flavoring_agents", "fragrances")
OUTPUT_ROOT = ROOT / "results" / "merged_flavor_fragrance"
RUN_ROOT = ROOT / "app" / "output" / "merged_flavor_fragrance"


def canonical_set(category: str) -> set[str]:
    output: set[str] = set()
    for smiles in category_positive_smiles(category):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None:
            output.add(Chem.MolToSmiles(molecule, canonical=True))
    return output


def split(values: set[str], fold: int) -> tuple[list[str], list[str]]:
    train, holdout = [], []
    for value in sorted(values):
        bucket = int(hashlib.sha256(value.encode()).hexdigest()[:8], 16) % 5
        (holdout if bucket == fold else train).append(value)
    return train, holdout


def write_smiles(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SMILES"])
        writer.writeheader()
        writer.writerows({"SMILES": value} for value in values)


def score(config: dict, values: list[str]) -> np.ndarray:
    rows = _choi_prepare_rows(values)
    descriptor, structural = _choi_two_scores(
        rows, config.get("selected_props", []), config.get("ranges", {}), config.get("pattern_weights", {})
    )
    weight = float(config.get("best_w", 0.5))
    return weight * descriptor + (1 - weight) * structural


def metrics(config: dict, positives: list[str], negatives: list[str]) -> dict[str, float]:
    positive = score(config, positives)
    negative = score(config, negatives)
    labels = np.r_[np.ones(len(positive)), np.zeros(len(negative))]
    values = np.r_[positive, negative]
    auc = float(roc_auc_score(labels, values))
    fpr, tpr, thresholds = roc_curve(labels, values)
    balanced = (tpr + 1 - fpr) / 2
    index = int(np.argmax(balanced))
    return {
        "auc": auc,
        "balanced_accuracy": float(balanced[index]),
        "threshold": float(thresholds[index]),
        "tpr": float(tpr[index]),
        "fpr": float(fpr[index]),
    }


def evaluate_fold(fold: int, bayes_trials: int) -> dict:
    targets = {category: canonical_set(category) for category in TARGET_CATEGORIES}
    target = set().union(*targets.values())
    sources = {category: canonical_set(category) - target for category in SOURCE_CATEGORIES}
    positive_train, positive_holdout = split(target, fold)
    negative_train, negative_holdout = [], []
    source_holdouts: dict[str, list[str]] = {}
    for category, values in sources.items():
        train, holdout = split(values, fold)
        negative_train.extend(train)
        negative_holdout.extend(holdout)
        source_holdouts[category] = holdout
    run_dir = RUN_ROOT / f"fold_{fold}"
    positive_path = run_dir / "positive.csv"
    negative_path = run_dir / "negative.csv"
    model_path = run_dir / "candidate.json"
    write_smiles(positive_path, positive_train)
    write_smiles(negative_path, sorted(set(negative_train)))
    build_choi_model(
        positive_csv=positive_path,
        negative_source_csv=negative_path,
        model_id="final_flavor_fragrance",
        label="Flavor and fragrance category score",
        category="flavor_fragrance",
        output_path=model_path,
        tanimoto_threshold=1.01,
        bayes_trials=bayes_trials,
        seed=20260816 + fold,
        category_prior=None,
    )
    config = json.loads(model_path.read_text(encoding="utf-8"))
    config["threshold"] = metrics(config, positive_train, sorted(set(negative_train)))["threshold"]
    overall = metrics(config, positive_holdout, sorted(set(negative_holdout)))
    source_response = {
        category: float(np.mean(score(config, values) >= config["threshold"]))
        for category, values in source_holdouts.items()
        if values
    }
    result = {
        "fold": fold,
        "positive_holdout_count": len(positive_holdout),
        "negative_holdout_count": len(set(negative_holdout)),
        "metrics": overall,
        "source_threshold_response": source_response,
        "worst_source": max(source_response, key=source_response.get),
        "worst_source_response": max(source_response.values()),
        "positive_response": float(np.mean(score(config, positive_holdout) >= config["threshold"])),
        "model_path": str(model_path),
    }
    out = OUTPUT_ROOT / f"fold_{fold}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", type=int, default=3)
    parser.add_argument("--bayes-trials", type=int, default=30)
    args = parser.parse_args()
    results = [evaluate_fold(fold, args.bayes_trials) for fold in range(args.folds)]
    targets = set().union(*(canonical_set(category) for category in TARGET_CATEGORIES))
    negatives = set().union(*(canonical_set(category) for category in SOURCE_CATEGORIES)) - targets
    positive_path = RUN_ROOT / "final_full" / "positive.csv"
    negative_path = RUN_ROOT / "final_full" / "negative.csv"
    model_path = RUN_ROOT / "final_full" / "candidate.json"
    write_smiles(positive_path, sorted(targets))
    write_smiles(negative_path, sorted(negatives))
    build_choi_model(
        positive_csv=positive_path,
        negative_source_csv=negative_path,
        model_id="final_flavor_fragrance",
        label="Flavor and fragrance category score",
        category="flavor_fragrance",
        output_path=model_path,
        tanimoto_threshold=1.01,
        bayes_trials=args.bayes_trials,
        seed=20260816,
        category_prior=None,
    )
    final_model = json.loads(model_path.read_text(encoding="utf-8"))
    final_model["threshold"] = metrics(final_model, sorted(targets), sorted(negatives))["threshold"]
    final_model["merged_categories"] = list(TARGET_CATEGORIES)
    final_model["validation_summary"] = results
    final_model["description"] = (
        "Merged flavor-and-fragrance score trained after the source categories showed 46% exact structural overlap "
        "and reciprocal cross-response; evaluated against exact-overlap-excluded cross-category negatives."
    )
    model_path.write_text(json.dumps(final_model, ensure_ascii=False, indent=2), encoding="utf-8")
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "summary.json").write_text(
        json.dumps({"folds": results, "final_model": str(model_path)}, indent=2), encoding="utf-8"
    )
    for result in results:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
