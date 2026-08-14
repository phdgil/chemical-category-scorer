from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from rdkit import Chem
from sklearn.metrics import roc_auc_score

from app.qed_inspired_validation import MODEL_IDS, ROOT, category_positive_smiles, score_vector

CANDIDATES = (
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


def canonical_smiles(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(molecule, canonical=True) if molecule is not None else None


def positive_sets() -> dict[str, set[str]]:
    output: dict[str, set[str]] = {}
    for category in CANDIDATES:
        values = {
            canonical
            for smiles in category_positive_smiles(category)
            if (canonical := canonical_smiles(smiles)) is not None
        }
        output[category] = values
    return output


def threshold_fraction(smiles: set[str], model_id: str) -> float:
    scores, _descriptor, _structural, threshold, _valid = score_vector(sorted(smiles), model_id)
    return float(np.mean(scores >= threshold))


def build_rows() -> list[dict[str, object]]:
    sets = positive_sets()
    rows: list[dict[str, object]] = []
    for category in CANDIDATES:
        model_id = MODEL_IDS.get(category, f"final_{category}")
        own = sets[category]
        own_response = threshold_fraction(own, model_id)
        source_responses: list[tuple[float, str, int]] = []
        for source in CANDIDATES:
            if source == category:
                continue
            source_exclusive = sets[source] - own
            source_responses.append(
                (threshold_fraction(source_exclusive, model_id), source, len(source_exclusive))
            )
        worst_response, worst_source, worst_count = max(source_responses)
        hard_background = set().union(*(sets[source] for source in CANDIDATES if source != category)) - own
        positive_scores, *_ = score_vector(sorted(own), model_id)
        negative_scores, *_ = score_vector(sorted(hard_background), model_id)
        labels = np.concatenate((np.ones(len(positive_scores)), np.zeros(len(negative_scores))))
        scores = np.concatenate((positive_scores, negative_scores))
        hard_auc = float(roc_auc_score(labels, scores))
        retained = own_response > worst_response and hard_auc > 0.5
        rows.append(
            {
                "category": category,
                "own_positive_count": len(own),
                "own_threshold_response": own_response,
                "worst_exact_overlap_excluded_source": worst_source,
                "worst_source_exclusive_count": worst_count,
                "worst_source_threshold_response": worst_response,
                "own_minus_worst_response": own_response - worst_response,
                "hard_background_count": len(hard_background),
                "exact_overlap_excluded_hard_auc": hard_auc,
                "publication_decision": "retain" if retained else "exclude",
                "decision_rule": "retain only if own response exceeds every exact-overlap-excluded source response and hard-background AUC exceeds 0.5",
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "qed_inspired_analysis" / "publication_specificity_revalidation.csv",
    )
    args = parser.parse_args()
    rows = build_rows()
    write_csv(args.output, rows)
    for row in rows:
        print(
            f"{row['category']}: {row['publication_decision']} "
            f"own={100 * float(row['own_threshold_response']):.1f}% "
            f"worst={100 * float(row['worst_source_threshold_response']):.1f}% "
            f"hard_auc={float(row['exact_overlap_excluded_hard_auc']):.3f}"
        )
    print(args.output)


if __name__ == "__main__":
    main()
