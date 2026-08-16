from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from sklearn.metrics import roc_auc_score, roc_curve

from app.algorithm_score_engine import _prepare_runtime_model, _score_han, score_smiles
from app.build_scoring_models import build_choi_model
from app.qed_inspired_validation import ROOT, category_positive_smiles, retained_negatives
from app.sequential_hard_negative_rebuild import canonical_set, metrics as choi_metrics, score_config, write_smiles

ORIGINAL_CATEGORIES = (
    "animal_drugs",
    "human_drugs",
    "cosmetics",
    "endocrine_disruptors",
    "flavoring_agents",
    "food_additives",
    "food_contact_substances",
    "fragrances",
    "pesticides",
    "solvents",
    "surfactants",
)
TARGETS = (
    "endocrine_disruptors",
    "flavor_fragrance",
    "pesticides",
    "surfactants",
    "animal_drugs",
    "human_drugs",
    "food_additives",
    "solvents",
)
EXTERNAL = {
    "endocrine_disruptors": [
        (ROOT / "results/external_validation/deduct_v3_endocrine/endocrine_disruptors_true_external_candidates_deduct_I-III.csv", "parent_smiles", None),
    ],
    "flavor_fragrance": [
        (ROOT / "results/external_validation/resolved/flavoring_agents_eu.csv", "SMILES", "resolved"),
    ],
    "pesticides": [
        (ROOT / "results/external_validation/analysis/pesticides/pesticides_true_external_candidates.csv", "parent_smiles", None),
    ],
    "surfactants": [
        (ROOT / "results/external_validation/resolved/surfactants_epa_scil.csv", "SMILES", "resolved"),
    ],
    "animal_drugs": [
        (ROOT / "results/external_validation/resolved/animal_drugs_canada_dpd.csv", "SMILES", "resolved"),
    ],
    "human_drugs": [
        (ROOT / "results/external_validation/analysis/human_drugs/human_drugs_true_external_candidates.csv", "parent_smiles", None),
        (ROOT / "results/external_validation/resolved/human_drugs_canada_dpd.csv", "SMILES", "resolved"),
    ],
    "food_additives": [
        (ROOT / "results/external_validation/resolved/food_additives_canada.csv", "SMILES", "resolved"),
    ],
    "solvents": [
        (ROOT / "results/external_validation/resolved/solvents_epa_scil.csv", "SMILES", "resolved"),
    ],
}
MODEL_IDS = {
    "endocrine_disruptors": "han_endocrine_disruptors",
    "flavor_fragrance": "final_flavor_fragrance",
    "pesticides": "final_pesticides",
    "surfactants": "final_surfactants",
    "animal_drugs": "final_animal_drugs",
    "human_drugs": "final_human_drugs",
    "food_additives": "final_food_additives",
    "solvents": "final_solvents",
}
OUTPUT_ROOT = ROOT / "results/combined_external_positive_rebuild"
RUN_ROOT = ROOT / "app/output/combined_external_positive_rebuild"


def read_external(category: str) -> set[str]:
    values: list[str] = []
    for path, column, required_status in EXTERNAL[category]:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                if required_status and row.get("resolution_status") != required_status:
                    continue
                values.append(row.get(column, ""))
    return canonical_set(values)


def pubchem_positive(category: str) -> set[str]:
    if category == "flavor_fragrance":
        return canonical_set(category_positive_smiles("flavoring_agents")) | canonical_set(
            category_positive_smiles("fragrances")
        )
    return canonical_set(category_positive_smiles(category))


def bucket(value: str) -> int:
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:8], 16) % 3


def threshold_metrics(positive: np.ndarray, negative: np.ndarray) -> dict[str, float]:
    labels = np.r_[np.ones(len(positive)), np.zeros(len(negative))]
    scores = np.r_[positive, negative]
    auc = float(roc_auc_score(labels, scores))
    fpr, tpr, thresholds = roc_curve(labels, scores)
    balanced = (tpr + 1 - fpr) / 2
    index = int(np.argmax(balanced))
    return {
        "auc": auc,
        "threshold": float(thresholds[index]),
        "balanced_accuracy": float(balanced[index]),
        "tpr": float(tpr[index]),
        "fpr": float(fpr[index]),
    }


def baseline_scores(category: str, values: list[str]) -> np.ndarray:
    return np.asarray([score_smiles(value, MODEL_IDS[category]).score for value in values], dtype=float)


def descriptor_values(values: list[str]) -> dict[str, np.ndarray]:
    output: dict[str, list[float]] = {name: [] for name in ("MW", "LogP", "HBA", "HBD", "TPSA", "RotB", "ArRings")}
    for value in values:
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            continue
        row = {
            "MW": Descriptors.MolWt(molecule),
            "LogP": Crippen.MolLogP(molecule),
            "HBA": rdMolDescriptors.CalcNumHBA(molecule),
            "HBD": rdMolDescriptors.CalcNumHBD(molecule),
            "TPSA": rdMolDescriptors.CalcTPSA(molecule),
            "RotB": rdMolDescriptors.CalcNumRotatableBonds(molecule),
            "ArRings": rdMolDescriptors.CalcNumAromaticRings(molecule),
        }
        for name, number in row.items():
            output[name].append(float(number))
    return {name: np.asarray(numbers) for name, numbers in output.items()}


def build_han_candidate(positives: list[str], output_dir: Path) -> tuple[dict[str, Any], Callable[[list[str]], np.ndarray]]:
    source = json.loads((ROOT / "app/data/models/han_endocrine_disruptors.json").read_text(encoding="utf-8"))
    reference_path = output_dir / "combined_reference.csv"
    with reference_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["canonical_smiles"])
        writer.writeheader()
        writer.writerows({"canonical_smiles": value} for value in positives)
    source["reference_smiles_path"] = str(reference_path.resolve())
    source["model_id"] = "candidate_combined_endocrine_disruptors"
    source["description"] = "Han endocrine score with descriptor statistics and similarity references updated from combined PubChem and external positives."
    source["descriptor_stats"] = {}
    for name, values in descriptor_values(positives).items():
        mean = float(np.mean(values))
        lower = values[values < mean]
        upper = values[values >= mean]
        source["descriptor_stats"][name] = {
            "mu": mean,
            "sig_down": float(np.mean(mean - lower)) if len(lower) else 1.0,
            "sig_up": float(np.mean(upper - mean)) if len(upper) else 1.0,
        }
    runtime = _prepare_runtime_model(source)

    def scorer(values: list[str]) -> np.ndarray:
        scores: list[float] = []
        for value in values:
            molecule = Chem.MolFromSmiles(value)
            if molecule is not None:
                scores.append(_score_han(value, molecule, runtime, source).score)
        return np.asarray(scores)

    return source, scorer


def build_negative_sets(target: str, target_positive: set[str], fold: int) -> tuple[list[str], list[str]]:
    train: set[str] = set()
    holdout: set[str] = set()
    for category in ORIGINAL_CATEGORIES:
        if category == target or (target == "flavor_fragrance" and category in {"flavoring_agents", "fragrances"}):
            continue
        values = canonical_set(category_positive_smiles(category))
        external_key = "flavor_fragrance" if category == "flavoring_agents" else category
        if external_key in EXTERNAL:
            values |= read_external(external_key)
        values -= target_positive
        for value in values:
            (holdout if bucket(value) == fold else train).add(value)
    return sorted(train), sorted(holdout)


def evaluate_target(category: str, fold: int, bayes_trials: int) -> dict[str, Any]:
    pubchem = pubchem_positive(category)
    external = read_external(category) - pubchem
    external_train = sorted(value for value in external if bucket(value) != fold)
    external_holdout = sorted(value for value in external if bucket(value) == fold)
    combined_train = sorted(pubchem | set(external_train))
    target_all = pubchem | external
    negative_train, negative_holdout = build_negative_sets(category, target_all, fold)
    output_dir = RUN_ROOT / category / f"fold_{fold}"
    output_dir.mkdir(parents=True, exist_ok=True)

    if category == "endocrine_disruptors":
        candidate, candidate_score = build_han_candidate(combined_train, output_dir)
        train_positive_scores = candidate_score(combined_train)
        train_negative_scores = candidate_score(negative_train)
        candidate["threshold"] = threshold_metrics(train_positive_scores, train_negative_scores)["threshold"]
    else:
        positive_path = output_dir / "combined_positive_train.csv"
        negative_path = output_dir / "combined_negative_train.csv"
        model_path = output_dir / "candidate.json"
        write_smiles(positive_path, combined_train)
        write_smiles(negative_path, negative_train)
        build_choi_model(
            positive_csv=positive_path,
            negative_source_csv=negative_path,
            model_id=f"candidate_combined_{category}",
            label=f"Combined-source {category.replace('_', ' ')} candidate",
            category=category,
            output_path=model_path,
            tanimoto_threshold=1.01,
            bayes_trials=bayes_trials,
            seed=20260816 + fold,
            category_prior=None,
        )
        candidate = json.loads(model_path.read_text(encoding="utf-8"))
        candidate_score = lambda values: score_config(candidate, values)

    candidate_external = candidate_score(external_holdout)
    baseline_external = baseline_scores(category, external_holdout)
    candidate_negative = candidate_score(negative_holdout)
    baseline_negative = baseline_scores(category, negative_holdout)
    candidate_external_metrics = threshold_metrics(candidate_external, candidate_negative)
    baseline_external_metrics = threshold_metrics(baseline_external, baseline_negative)

    retained = retained_negatives(category, ROOT / "results/qed_inspired_analysis")
    candidate_original = threshold_metrics(candidate_score(sorted(pubchem)), candidate_score(retained))
    baseline_original = threshold_metrics(baseline_scores(category, sorted(pubchem)), baseline_scores(category, retained))
    result = {
        "category": category,
        "fold": fold,
        "pubchem_positive_count": len(pubchem),
        "external_total_count": len(external),
        "external_train_count": len(external_train),
        "external_holdout_count": len(external_holdout),
        "negative_train_count": len(negative_train),
        "negative_holdout_count": len(negative_holdout),
        "baseline_external_holdout": baseline_external_metrics,
        "candidate_external_holdout": candidate_external_metrics,
        "external_holdout_auc_delta": candidate_external_metrics["auc"] - baseline_external_metrics["auc"],
        "baseline_original_benchmark": baseline_original,
        "candidate_original_benchmark": candidate_original,
        "original_auc_delta": candidate_original["auc"] - baseline_original["auc"],
        "candidate_config": candidate,
    }
    result_path = OUTPUT_ROOT / category / f"fold_{fold}.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", default=",".join(TARGETS))
    parser.add_argument("--folds", type=int, default=3)
    parser.add_argument("--bayes-trials", type=int, default=20)
    args = parser.parse_args()
    targets = [value.strip() for value in args.targets.split(",") if value.strip()]
    results = [evaluate_target(target, fold, args.bayes_trials) for target in targets for fold in range(args.folds)]
    summary: dict[str, Any] = {}
    for target in targets:
        selected = [row for row in results if row["category"] == target]
        summary[target] = {
            "folds": selected,
            "mean_external_holdout_auc_delta": float(np.mean([row["external_holdout_auc_delta"] for row in selected])),
            "minimum_external_holdout_auc_delta": float(np.min([row["external_holdout_auc_delta"] for row in selected])),
            "mean_original_auc_delta": float(np.mean([row["original_auc_delta"] for row in selected])),
            "promotion_gate": bool(
                min(row["external_holdout_auc_delta"] for row in selected) > 0
                and np.mean([row["external_holdout_auc_delta"] for row in selected]) >= 0.02
                and np.mean([row["original_auc_delta"] for row in selected]) >= -0.02
            ),
        }
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    for target, values in summary.items():
        print(target, {key: value for key, value in values.items() if key != "folds"})


if __name__ == "__main__":
    main()
