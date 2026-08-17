from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy import sparse
from sklearn.metrics import roc_auc_score, roc_curve

from app.qed_inspired_validation import ROOT, category_positive_smiles, retained_negatives, score_vector
from app.sequential_hard_negative_rebuild import canonical_set

OUTPUT_ROOT = ROOT / "results/network_random_walk_pattern_rebuild"
TARGETS = ("endocrine_disruptors", "flavor_fragrance", "pesticides", "surfactants")
MODEL_IDS = {
    "endocrine_disruptors": "han_endocrine_disruptors",
    "flavor_fragrance": "final_flavor_fragrance",
    "pesticides": "final_pesticides",
    "surfactants": "final_surfactants",
}
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


def target_positives(target: str) -> set[str]:
    if target == "flavor_fragrance":
        return canonical_set(category_positive_smiles("flavoring_agents")) | canonical_set(
            category_positive_smiles("fragrances")
        )
    return canonical_set(category_positive_smiles(target))


def hard_negatives(target: str, positives: set[str]) -> set[str]:
    excluded = {target}
    if target == "flavor_fragrance":
        excluded = {"flavoring_agents", "fragrances"}
    values: set[str] = set()
    for category in ORIGINAL_CATEGORIES:
        if category not in excluded:
            values |= canonical_set(category_positive_smiles(category))
    return values - positives


def fold_bucket(smiles: str) -> int:
    return int(hashlib.sha256(smiles.encode("utf-8")).hexdigest()[:8], 16) % 3


def pattern_scores(
    values: list[str],
    patterns: list[dict[str, float | int | str]],
    cache: dict[str, set[str]],
) -> np.ndarray:
    total = sum(float(pattern["weight"]) for pattern in patterns)
    if total <= 0:
        return np.zeros(len(values), dtype=float)
    weights = {str(pattern["fragment_smiles"]): float(pattern["weight"]) for pattern in patterns}
    return np.asarray(
        [sum(weights.get(fragment, 0.0) for fragment in cache[value]) / total for value in values],
        dtype=float,
    )


def auc(positive: np.ndarray, negative: np.ndarray) -> float:
    labels = np.r_[np.ones(len(positive)), np.zeros(len(negative))]
    return float(roc_auc_score(labels, np.r_[positive, negative]))


def balanced_threshold(positive: np.ndarray, negative: np.ndarray) -> tuple[float, float]:
    labels = np.r_[np.ones(len(positive)), np.zeros(len(negative))]
    scores = np.r_[positive, negative]
    fpr, tpr, thresholds = roc_curve(labels, scores)
    balanced = (tpr + 1 - fpr) / 2
    index = int(np.argmax(balanced))
    return float(thresholds[index]), float(balanced[index])


def optimize_mix(
    baseline_positive: np.ndarray,
    baseline_negative: np.ndarray,
    pattern_positive: np.ndarray,
    pattern_negative: np.ndarray,
) -> tuple[float, float, float]:
    best: tuple[float, float, float, float] | None = None
    for baseline_weight in np.linspace(0, 1, 21):
        positive = baseline_weight * baseline_positive + (1 - baseline_weight) * pattern_positive
        negative = baseline_weight * baseline_negative + (1 - baseline_weight) * pattern_negative
        candidate_auc = auc(positive, negative)
        threshold, balanced = balanced_threshold(positive, negative)
        key = (candidate_auc, balanced, baseline_weight, threshold)
        if best is None or key[:2] > best[:2]:
            best = key
    assert best is not None
    return best[2], best[3], best[0]


def neighborhood_fragments(smiles: str, minimum_atoms: int, maximum_radius: int) -> set[str]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return set()
    output: set[str] = set()
    for atom in molecule.GetAtoms():
        center = atom.GetIdx()
        for radius in range(1, maximum_radius + 1):
            bonds = list(Chem.FindAtomEnvironmentOfRadiusN(molecule, radius, center))
            if not bonds:
                continue
            atoms = {center}
            for bond_index in bonds:
                bond = molecule.GetBondWithIdx(bond_index)
                atoms.add(bond.GetBeginAtomIdx())
                atoms.add(bond.GetEndAtomIdx())
            try:
                fragment = Chem.MolFragmentToSmiles(
                    molecule,
                    atomsToUse=sorted(atoms),
                    bondsToUse=bonds,
                    canonical=True,
                    isomericSmiles=False,
                )
            except RuntimeError:
                continue
            fragment_molecule = Chem.MolFromSmiles(fragment)
            if fragment_molecule is not None and fragment_molecule.GetNumHeavyAtoms() >= minimum_atoms:
                output.add(fragment)
    return output


def build_cache(values: set[str], minimum_atoms: int, maximum_radius: int) -> dict[str, set[str]]:
    return {
        value: neighborhood_fragments(value, minimum_atoms, maximum_radius)
        for value in values
    }


def bipartite_pagerank(
    incidence: sparse.csr_matrix,
    seed_rows: np.ndarray,
    damping: float,
    tolerance: float,
    maximum_iterations: int,
) -> np.ndarray:
    molecule_count, fragment_count = incidence.shape
    molecule_degree = np.asarray(incidence.sum(axis=1)).ravel()
    fragment_degree = np.asarray(incidence.sum(axis=0)).ravel()
    molecule_degree[molecule_degree == 0] = 1
    fragment_degree[fragment_degree == 0] = 1
    restart = np.zeros(molecule_count, dtype=float)
    restart[seed_rows] = 1.0 / len(seed_rows)
    molecule_probability = restart.copy()
    fragment_probability = np.zeros(fragment_count, dtype=float)
    for _ in range(maximum_iterations):
        next_fragment = damping * np.asarray(
            incidence.T @ (molecule_probability / molecule_degree)
        ).ravel()
        next_molecule = (1 - damping) * restart + damping * np.asarray(
            incidence @ (fragment_probability / fragment_degree)
        ).ravel()
        difference = np.abs(next_molecule - molecule_probability).sum() + np.abs(
            next_fragment - fragment_probability
        ).sum()
        molecule_probability = next_molecule
        fragment_probability = next_fragment
        if difference < tolerance:
            break
    total = fragment_probability.sum()
    return fragment_probability / total if total > 0 else fragment_probability


def discover_network_patterns(
    positive_values: list[str],
    negative_values: list[str],
    cache: dict[str, set[str]],
    damping: float,
    minimum_prevalence: float,
    minimum_enrichment: float,
    maximum_patterns: int,
) -> list[dict[str, float | int | str]]:
    molecules = positive_values + negative_values
    positive_counts: Counter[str] = Counter()
    negative_counts: Counter[str] = Counter()
    for value in positive_values:
        positive_counts.update(cache[value])
    for value in negative_values:
        negative_counts.update(cache[value])
    minimum_count = max(10, math.ceil(minimum_prevalence * len(positive_values)))
    vocabulary = sorted(
        fragment
        for fragment, count in positive_counts.items()
        if count >= minimum_count
    )
    fragment_index = {fragment: index for index, fragment in enumerate(vocabulary)}
    row_indices: list[int] = []
    column_indices: list[int] = []
    for row_index, value in enumerate(molecules):
        for fragment in cache[value]:
            column_index = fragment_index.get(fragment)
            if column_index is not None:
                row_indices.append(row_index)
                column_indices.append(column_index)
    incidence = sparse.csr_matrix(
        (np.ones(len(row_indices), dtype=float), (row_indices, column_indices)),
        shape=(len(molecules), len(vocabulary)),
    )
    positive_walk = bipartite_pagerank(
        incidence,
        np.arange(len(positive_values)),
        damping,
        1e-10,
        200,
    )
    negative_walk = bipartite_pagerank(
        incidence,
        np.arange(len(positive_values), len(molecules)),
        damping,
        1e-10,
        200,
    )
    epsilon = 1e-15
    rows: list[dict[str, float | int | str]] = []
    for fragment, index in fragment_index.items():
        positive_count = positive_counts[fragment]
        negative_count = negative_counts.get(fragment, 0)
        positive_prevalence = positive_count / len(positive_values)
        negative_prevalence = negative_count / len(negative_values)
        prevalence_enrichment = ((positive_count + 1) / (len(positive_values) + 2)) / (
            (negative_count + 1) / (len(negative_values) + 2)
        )
        walk_enrichment = (positive_walk[index] + epsilon) / (negative_walk[index] + epsilon)
        if walk_enrichment < minimum_enrichment or prevalence_enrichment < 1.25:
            continue
        weight = math.log2(walk_enrichment) * math.sqrt(positive_prevalence)
        rows.append(
            {
                "fragment_smiles": fragment,
                "positive_count": positive_count,
                "negative_count": negative_count,
                "positive_prevalence": positive_prevalence,
                "negative_prevalence": negative_prevalence,
                "prevalence_enrichment": prevalence_enrichment,
                "positive_walk_probability": float(positive_walk[index]),
                "negative_walk_probability": float(negative_walk[index]),
                "walk_enrichment": float(walk_enrichment),
                "weight": weight,
            }
        )
    rows.sort(key=lambda row: (float(row["weight"]), float(row["positive_prevalence"])), reverse=True)
    return rows[:maximum_patterns]


def evaluate_target(
    target: str,
    minimum_atoms: int,
    maximum_radius: int,
    damping: float,
    minimum_prevalence: float,
    minimum_enrichment: float,
    maximum_patterns: int,
) -> list[dict[str, object]]:
    positives = target_positives(target)
    negatives = hard_negatives(target, positives)
    retained = canonical_set(retained_negatives(target, ROOT / "results/qed_inspired_analysis"))
    all_values = positives | negatives | retained
    cache = build_cache(all_values, minimum_atoms, maximum_radius)
    ordered_values = sorted(all_values)
    combined, _descriptor, _structural, _threshold, _valid = score_vector(ordered_values, MODEL_IDS[target])
    baseline_lookup = dict(zip(ordered_values, combined))

    def baseline(values: list[str]) -> np.ndarray:
        return np.asarray([baseline_lookup[value] for value in values], dtype=float)

    results: list[dict[str, object]] = []
    for fold in range(3):
        positive_train = sorted(value for value in positives if fold_bucket(value) != fold)
        positive_holdout = sorted(value for value in positives if fold_bucket(value) == fold)
        negative_train = sorted(value for value in negatives if fold_bucket(value) != fold)
        negative_holdout = sorted(value for value in negatives if fold_bucket(value) == fold)
        retained_holdout = sorted(value for value in retained if fold_bucket(value) == fold)
        patterns = discover_network_patterns(
            positive_train,
            negative_train,
            cache,
            damping,
            minimum_prevalence,
            minimum_enrichment,
            maximum_patterns,
        )
        train_baseline_positive = baseline(positive_train)
        train_baseline_negative = baseline(negative_train)
        train_pattern_positive = pattern_scores(positive_train, patterns, cache)
        train_pattern_negative = pattern_scores(negative_train, patterns, cache)
        baseline_weight, threshold, training_auc = optimize_mix(
            train_baseline_positive,
            train_baseline_negative,
            train_pattern_positive,
            train_pattern_negative,
        )
        holdout_baseline_positive = baseline(positive_holdout)
        holdout_baseline_negative = baseline(negative_holdout)
        holdout_pattern_positive = pattern_scores(positive_holdout, patterns, cache)
        holdout_pattern_negative = pattern_scores(negative_holdout, patterns, cache)
        holdout_candidate_positive = baseline_weight * holdout_baseline_positive + (1 - baseline_weight) * holdout_pattern_positive
        holdout_candidate_negative = baseline_weight * holdout_baseline_negative + (1 - baseline_weight) * holdout_pattern_negative
        retained_baseline_negative = baseline(retained_holdout)
        retained_pattern_negative = pattern_scores(retained_holdout, patterns, cache)
        retained_candidate_negative = baseline_weight * retained_baseline_negative + (1 - baseline_weight) * retained_pattern_negative
        baseline_hard_auc = auc(holdout_baseline_positive, holdout_baseline_negative)
        candidate_hard_auc = auc(holdout_candidate_positive, holdout_candidate_negative)
        baseline_original_auc = auc(holdout_baseline_positive, retained_baseline_negative)
        candidate_original_auc = auc(holdout_candidate_positive, retained_candidate_negative)
        results.append(
            {
                "category": target,
                "fold": fold,
                "selected_pattern_count": len(patterns),
                "baseline_weight": baseline_weight,
                "network_pattern_weight": 1 - baseline_weight,
                "candidate_threshold": threshold,
                "training_auc": training_auc,
                "baseline_hard_auc": baseline_hard_auc,
                "candidate_hard_auc": candidate_hard_auc,
                "hard_auc_delta": candidate_hard_auc - baseline_hard_auc,
                "baseline_original_auc": baseline_original_auc,
                "candidate_original_auc": candidate_original_auc,
                "original_auc_delta": candidate_original_auc - baseline_original_auc,
                "patterns": patterns,
            }
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", default=",".join(TARGETS))
    parser.add_argument("--minimum-atoms", type=int, default=3)
    parser.add_argument("--maximum-radius", type=int, default=3)
    parser.add_argument("--damping", type=float, default=0.85)
    parser.add_argument("--minimum-prevalence", type=float, default=0.01)
    parser.add_argument("--minimum-enrichment", type=float, default=1.5)
    parser.add_argument("--maximum-patterns", type=int, default=24)
    args = parser.parse_args()
    summary: dict[str, object] = {}
    for target in [item.strip() for item in args.targets.split(",") if item.strip()]:
        folds = evaluate_target(
            target,
            args.minimum_atoms,
            args.maximum_radius,
            args.damping,
            args.minimum_prevalence,
            args.minimum_enrichment,
            args.maximum_patterns,
        )
        hard_deltas = [float(row["hard_auc_delta"]) for row in folds]
        original_deltas = [float(row["original_auc_delta"]) for row in folds]
        promotion = min(hard_deltas) > 0 and np.mean(hard_deltas) >= 0.02 and np.mean(original_deltas) >= -0.02
        summary[target] = {
            "folds": folds,
            "mean_hard_auc_delta": float(np.mean(hard_deltas)),
            "minimum_hard_auc_delta": float(np.min(hard_deltas)),
            "mean_original_auc_delta": float(np.mean(original_deltas)),
            "promotion_gate": bool(promotion),
        }
        print(target, {key: value for key, value in summary[target].items() if key != "folds"})
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
