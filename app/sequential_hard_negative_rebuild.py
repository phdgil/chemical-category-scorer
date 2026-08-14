from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from sklearn.metrics import roc_auc_score, roc_curve

from app.build_scoring_models import _choi_prepare_rows, _choi_two_scores, build_choi_model
from app.qed_inspired_validation import CATEGORIES, MODEL_IDS, ROOT, category_positive_smiles, read_csv, write_csv

TARGETS = (
    "animal_drugs",
    "human_drugs",
    "food_contact_substances",
    "cosmetics",
    "food_additives",
)
OUTPUT_ROOT = ROOT / "results" / "sequential_hard_negative_rebuild"
RUN_ROOT = ROOT / "app" / "output" / "sequential_hard_negative_rebuild"
MODEL_DIR = ROOT / "app" / "data" / "models"
PROBE_SMILES = {
    "caffeine": "CN1C(=O)N(C)c2ncn(C)c2C1=O",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "ddt": "Clc1ccc(C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl)cc1",
    "bisphenol_a": "CC(c1ccc(O)cc1)(c1ccc(O)cc1)C",
    "vanillin": "COc1cc(C=O)ccc1O",
    "sds": "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]",
    "ethanol": "CCO",
}


def canonicalize(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    return Chem.MolToSmiles(molecule, canonical=True)


def canonical_set(values: list[str]) -> set[str]:
    return {canonical for value in values if (canonical := canonicalize(value)) is not None}


def split(values: set[str], holdout_bucket: int = 0) -> tuple[list[str], list[str]]:
    train: list[str] = []
    holdout: list[str] = []
    for value in sorted(values):
        bucket = int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:8], 16) % 5
        (holdout if bucket == holdout_bucket else train).append(value)
    return train, holdout


def sample(values: list[str], limit: int, salt: str) -> list[str]:
    if limit <= 0 or len(values) <= limit:
        return values
    ranked = sorted(values, key=lambda value: hashlib.sha256(f"{salt}|{value}".encode("utf-8")).hexdigest())
    return ranked[:limit]


def write_smiles(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SMILES"])
        writer.writeheader()
        writer.writerows({"SMILES": value} for value in values)


def score_config(config: dict[str, Any], smiles_values: list[str]) -> np.ndarray:
    rows = _choi_prepare_rows(smiles_values)
    property_scores, structure_scores = _choi_two_scores(
        rows,
        config.get("selected_props", []),
        config.get("ranges", {}),
        config.get("pattern_weights", {}),
    )
    weight = float(config.get("best_w", 0.5))
    return weight * property_scores + (1.0 - weight) * structure_scores


def metrics(config: dict[str, Any], positives: list[str], negatives: list[str]) -> dict[str, float]:
    positive_scores = score_config(config, positives)
    negative_scores = score_config(config, negatives)
    labels = np.r_[np.ones(len(positive_scores)), np.zeros(len(negative_scores))]
    scores = np.r_[positive_scores, negative_scores]
    auc = float(roc_auc_score(labels, scores))
    fpr, tpr, thresholds = roc_curve(labels, scores)
    balanced = 0.5 * (tpr + 1.0 - fpr)
    index = int(np.argmax(balanced))
    threshold = float(thresholds[index])
    valid = np.flatnonzero((fpr <= 0.10) & np.isfinite(thresholds))
    fpr10_index = int(valid[np.argmax(tpr[valid])]) if valid.size else 0
    deployed_threshold = float(config["threshold"])
    deployed_tpr = float(np.mean(positive_scores >= deployed_threshold))
    deployed_fpr = float(np.mean(negative_scores >= deployed_threshold))
    return {
        "auc": auc,
        "balanced_accuracy": float(balanced[index]),
        "balanced_threshold": threshold,
        "balanced_tpr": float(tpr[index]),
        "balanced_fpr": float(fpr[index]),
        "fpr10_threshold": float(thresholds[fpr10_index]),
        "fpr10_tpr": float(tpr[fpr10_index]),
        "fpr10_fpr": float(fpr[fpr10_index]),
        "deployed_threshold": deployed_threshold,
        "deployed_tpr": deployed_tpr,
        "deployed_fpr": deployed_fpr,
        "deployed_balanced_accuracy": 0.5 * (deployed_tpr + 1.0 - deployed_fpr),
        "positive_count": len(positive_scores),
        "negative_count": len(negative_scores),
    }


def build_hard_sets(
    target: str,
    max_train_per_category: int,
    max_holdout_per_category: int,
    holdout_bucket: int,
) -> dict[str, Any]:
    category_sets = {category: canonical_set(category_positive_smiles(category)) for category in CATEGORIES}
    target_set = category_sets[target]
    positive_train, positive_holdout = split(target_set, holdout_bucket)
    negative_train: list[str] = []
    negative_holdout: list[str] = []
    source_counts: dict[str, dict[str, int]] = {}
    for category in CATEGORIES:
        if category == target:
            continue
        eligible = category_sets[category] - target_set
        train, holdout = split(eligible, holdout_bucket)
        train = sample(train, max_train_per_category, f"{target}|{category}|train")
        holdout = sample(holdout, max_holdout_per_category, f"{target}|{category}|holdout")
        negative_train.extend(train)
        negative_holdout.extend(holdout)
        source_counts[category] = {"train": len(train), "holdout": len(holdout)}
    return {
        "positive_train": sorted(set(positive_train)),
        "positive_holdout": sorted(set(positive_holdout)),
        "negative_train": sorted(set(negative_train)),
        "negative_holdout": sorted(set(negative_holdout)),
        "source_counts": source_counts,
    }


def run_target(
    target: str,
    bayes_trials: int,
    max_train_per_category: int,
    max_holdout_per_category: int,
    holdout_bucket: int,
    seed: int,
) -> None:
    if target not in TARGETS:
        raise ValueError(f"Unsupported target {target!r}; choose from {TARGETS}")
    sets = build_hard_sets(target, max_train_per_category, max_holdout_per_category, holdout_bucket)
    target_dir = RUN_ROOT / target / f"fold_{holdout_bucket}"
    positive_csv = target_dir / "positive_train.csv"
    negative_csv = target_dir / "hard_cross_category_train.csv"
    candidate_path = target_dir / f"candidate_{target}.json"
    write_smiles(positive_csv, sets["positive_train"])
    write_smiles(negative_csv, sets["negative_train"])

    build_choi_model(
        positive_csv=positive_csv,
        negative_source_csv=negative_csv,
        model_id=f"candidate_hard_{target}",
        label=f"Hard-negative rebuilt {target.replace('_', ' ')}",
        category=target,
        output_path=candidate_path,
        tanimoto_threshold=1.01,
        bayes_trials=bayes_trials,
        seed=seed,
        category_prior="auto" if target in {"cosmetics", "food_contact_substances"} else None,
    )
    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    baseline = json.loads((MODEL_DIR / f"{MODEL_IDS[target]}.json").read_text(encoding="utf-8"))
    baseline_hard = metrics(baseline, sets["positive_holdout"], sets["negative_holdout"])
    candidate_hard = metrics(candidate, sets["positive_holdout"], sets["negative_holdout"])

    retained_negative_path = ROOT / "results" / "qed_inspired_analysis" / "retained_negatives" / f"{target}.csv"
    retained_negatives = [row["SMILES"] for row in read_csv(retained_negative_path)]
    full_positives = sorted(canonical_set(category_positive_smiles(target)))
    baseline_original = metrics(baseline, full_positives, retained_negatives)
    candidate_original = metrics(candidate, full_positives, retained_negatives)

    candidate["hard_cross_category_evaluation"] = candidate_hard
    candidate["baseline_hard_cross_category_evaluation"] = baseline_hard
    candidate["original_benchmark_evaluation"] = candidate_original
    candidate["baseline_original_benchmark_evaluation"] = baseline_original
    candidate["training_policy"] = "all other category structures excluding canonical target overlaps; no near-positive Tanimoto removal"
    candidate_path.write_text(json.dumps(candidate, ensure_ascii=False, indent=2), encoding="utf-8")

    probe_rows: list[dict[str, object]] = []
    for probe, smiles in PROBE_SMILES.items():
        canonical = canonicalize(smiles)
        if canonical is None:
            continue
        baseline_score = float(score_config(baseline, [canonical])[0])
        candidate_score = float(score_config(candidate, [canonical])[0])
        probe_rows.append(
            {
                "target": target,
                "probe": probe,
                "smiles": canonical,
                "baseline_score": baseline_score,
                "baseline_threshold": float(baseline["threshold"]),
                "baseline_likely": int(baseline_score >= float(baseline["threshold"])),
                "candidate_score": candidate_score,
                "candidate_threshold": float(candidate["threshold"]),
                "candidate_likely": int(candidate_score >= float(candidate["threshold"])),
            }
        )

    hard_auc_delta = candidate_hard["auc"] - baseline_hard["auc"]
    hard_ba_delta = candidate_hard["balanced_accuracy"] - baseline_hard["balanced_accuracy"]
    original_auc_delta = candidate_original["auc"] - baseline_original["auc"]
    promote = (
        hard_auc_delta >= 0.03
        and candidate_hard["deployed_balanced_accuracy"] >= baseline_hard["deployed_balanced_accuracy"] + 0.01
        and candidate_hard["deployed_fpr"] <= baseline_hard["deployed_fpr"] - 0.03
        and original_auc_delta >= -0.03
    )
    summary = {
        "target": target,
        "candidate_path": str(candidate_path),
        "positive_train_count": len(sets["positive_train"]),
        "positive_holdout_count": len(sets["positive_holdout"]),
        "hard_negative_train_count": len(sets["negative_train"]),
        "hard_negative_holdout_count": len(sets["negative_holdout"]),
        "baseline_hard": baseline_hard,
        "candidate_hard": candidate_hard,
        "hard_auc_delta": hard_auc_delta,
        "hard_balanced_accuracy_delta": hard_ba_delta,
        "baseline_original": baseline_original,
        "candidate_original": candidate_original,
        "original_auc_delta": original_auc_delta,
        "promotion_criteria_met": promote,
        "source_counts": sets["source_counts"],
    }
    output_dir = OUTPUT_ROOT / target / f"fold_{holdout_bucket}"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csv(output_dir / "probe_comparison.csv", probe_rows)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", choices=TARGETS, required=True)
    parser.add_argument("--bayes-trials", type=int, default=30)
    parser.add_argument("--max-train-per-category", type=int, default=800)
    parser.add_argument("--max-holdout-per-category", type=int, default=400)
    parser.add_argument("--holdout-bucket", type=int, choices=range(5), default=0)
    parser.add_argument("--seed", type=int, default=20260813)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_target(
        args.target,
        args.bayes_trials,
        args.max_train_per_category,
        args.max_holdout_per_category,
        args.holdout_bucket,
        args.seed,
    )


if __name__ == "__main__":
    main()
