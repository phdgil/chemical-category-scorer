from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve

from app.qed_inspired_validation import (
    CATEGORIES,
    DISPLAY,
    MODEL_IDS,
    ROOT,
    category_positive_smiles,
    score_vector,
    write_csv,
)

OUTPUT_DIR = ROOT / "results" / "qed_inspired_analysis"

EXAMPLE_SMILES = {
    "Caffeine": "CN1C(=O)N(C)c2ncn(C)c2C1=O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "DDT": "Clc1ccc(C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl)cc1",
    "Bisphenol A": "CC(c1ccc(O)cc1)(c1ccc(O)cc1)C",
    "Vanillin": "COc1cc(C=O)ccc1O",
    "SDS": "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]",
    "Ethanol": "CCO",
}


def threshold_at_fpr(labels: np.ndarray, scores: np.ndarray, maximum_fpr: float) -> float:
    fpr, tpr, thresholds = roc_curve(labels, scores)
    valid = np.flatnonzero((fpr <= maximum_fpr) & np.isfinite(thresholds))
    if not valid.size:
        return float(np.nextafter(np.max(scores), np.inf))
    best_tpr = np.max(tpr[valid])
    candidates = valid[tpr[valid] == best_tpr]
    return float(np.min(thresholds[candidates]))


def rates(labels: np.ndarray, scores: np.ndarray, threshold: float) -> tuple[float, float, float]:
    positive = labels == 1
    negative = ~positive
    tpr = float(np.mean(scores[positive] >= threshold))
    fpr = float(np.mean(scores[negative] >= threshold))
    balanced_accuracy = 0.5 * (tpr + 1.0 - fpr)
    return tpr, fpr, balanced_accuracy


def main() -> None:
    positive_sets = {category: set(category_positive_smiles(category)) for category in CATEGORIES}
    union_smiles = sorted(set().union(*positive_sets.values()))
    score_maps: dict[str, dict[str, float]] = {}
    original_thresholds: dict[str, float] = {}
    calibrated_thresholds: dict[str, float] = {}
    summary_rows: list[dict[str, object]] = []

    for category in CATEGORIES:
        values, _descriptor, _structural, original_threshold, valid = score_vector(
            union_smiles, MODEL_IDS[category]
        )
        score_map = dict(zip(valid, values, strict=True))
        score_maps[category] = score_map
        original_thresholds[category] = original_threshold
        usable = [smiles for smiles in union_smiles if smiles in score_map]
        labels = np.asarray([smiles in positive_sets[category] for smiles in usable], dtype=int)
        scores = np.asarray([score_map[smiles] for smiles in usable], dtype=float)
        calibrated_threshold = threshold_at_fpr(labels, scores, maximum_fpr=0.10)
        calibrated_thresholds[category] = calibrated_threshold
        old_tpr, old_fpr, old_ba = rates(labels, scores, original_threshold)
        new_tpr, new_fpr, new_ba = rates(labels, scores, calibrated_threshold)
        summary_rows.append(
            {
                "category": category,
                "model_id": MODEL_IDS[category],
                "union_count": len(usable),
                "positive_count": int(np.sum(labels)),
                "hard_cross_category_count": int(np.sum(labels == 0)),
                "auc_against_hard_cross_category": float(roc_auc_score(labels, scores)),
                "original_threshold": original_threshold,
                "original_tpr": old_tpr,
                "original_fpr": old_fpr,
                "original_balanced_accuracy": old_ba,
                "fpr10_threshold": calibrated_threshold,
                "fpr10_tpr": new_tpr,
                "fpr10_fpr": new_fpr,
                "fpr10_balanced_accuracy": new_ba,
            }
        )

    matrix_rows: list[dict[str, object]] = []
    for source_category in CATEGORIES:
        source = positive_sets[source_category]
        for score_category in CATEGORIES:
            score_map = score_maps[score_category]
            values = np.asarray([score_map[smiles] for smiles in source if smiles in score_map])
            threshold = calibrated_thresholds[score_category]
            matrix_rows.append(
                {
                    "source_category": source_category,
                    "score_category": score_category,
                    "positive_count": len(values),
                    "fpr10_threshold": threshold,
                    "fraction_at_or_above_fpr10_threshold": float(np.mean(values >= threshold)),
                }
            )

    example_rows: list[dict[str, object]] = []
    for compound, smiles in EXAMPLE_SMILES.items():
        for category in CATEGORIES:
            values, _descriptor, _structural, original_threshold, valid = score_vector(
                [smiles], MODEL_IDS[category]
            )
            if not valid:
                continue
            score = float(values[0])
            score_map = score_maps[category]
            positive_distribution = np.asarray(
                [score_map[value] for value in positive_sets[category] if value in score_map],
                dtype=float,
            )
            negative_distribution = np.asarray(
                [score_map[value] for value in union_smiles if value not in positive_sets[category] and value in score_map],
                dtype=float,
            )
            positive_tail = (float(np.sum(positive_distribution >= score)) + 0.5) / (
                len(positive_distribution) + 1.0
            )
            negative_tail = (float(np.sum(negative_distribution >= score)) + 0.5) / (
                len(negative_distribution) + 1.0
            )
            example_rows.append(
                {
                    "compound": compound,
                    "smiles": smiles,
                    "category": category,
                    "score": score,
                    "original_threshold": original_threshold,
                    "original_likely": int(score >= original_threshold),
                    "fpr10_threshold": calibrated_thresholds[category],
                    "fpr10_likely": int(score >= calibrated_thresholds[category]),
                    "positive_tail_fraction_at_query_score": positive_tail,
                    "hard_negative_tail_fraction_at_query_score": negative_tail,
                    "tail_likelihood_ratio": positive_tail / negative_tail,
                }
            )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT_DIR / "hard_cross_category_recalibration.csv", summary_rows)
    write_csv(OUTPUT_DIR / "hard_cross_category_fpr10_matrix.csv", matrix_rows)
    write_csv(OUTPUT_DIR / "coauthor_examples_recalibrated.csv", example_rows)
    calibration = {
        "schema_version": 2,
        "release_version": "2.1.0",
        "generated": "2026-08-17",
        "interpretation": {
            "shared": "score meets the original category threshold but not the cross-category-specific threshold",
            "high_specificity": "score meets the threshold calibrated to at most 10% response in hard cross-category structures",
            "below": "score is below the original category threshold",
        },
        "models": {
            str(row["model_id"]): {
                "category": str(row["category"]),
                "original_threshold": float(row["original_threshold"]),
                "high_specificity_threshold": float(row["fpr10_threshold"]),
                "hard_cross_category_auc": float(row["auc_against_hard_cross_category"]),
                "original_tpr": float(row["original_tpr"]),
                "original_hard_fpr": float(row["original_fpr"]),
                "high_specificity_tpr": float(row["fpr10_tpr"]),
                "high_specificity_fpr": float(row["fpr10_fpr"]),
            }
            for row in summary_rows
        },
    }
    calibration_path = ROOT / "app" / "data" / "cross_category_calibration.json"
    calibration_path.write_text(json.dumps(calibration, indent=2), encoding="utf-8")

    display = {category: DISPLAY[category].replace("\n", " ") for category in CATEGORIES}
    print(f"Union structures: {len(union_smiles):,}")
    for row in summary_rows:
        print(
            f"{display[str(row['category'])]}: FPR {100 * float(row['original_fpr']):.1f}% -> "
            f"{100 * float(row['fpr10_fpr']):.1f}%; TPR {100 * float(row['original_tpr']):.1f}% -> "
            f"{100 * float(row['fpr10_tpr']):.1f}%"
        )
    for compound in EXAMPLE_SMILES:
        rows = [row for row in example_rows if row["compound"] == compound]
        old = [display[str(row["category"])] for row in rows if row["original_likely"]]
        new = [display[str(row["category"])] for row in rows if row["fpr10_likely"]]
        ranked = sorted(rows, key=lambda row: float(row["tail_likelihood_ratio"]), reverse=True)
        top = [
            f"{display[str(row['category'])]} ({float(row['tail_likelihood_ratio']):.1f})"
            for row in ranked[:3]
        ]
        print(f"{compound}: original={old}; FPR10={new}; tail-LR top3={top}")


if __name__ == "__main__":
    main()
