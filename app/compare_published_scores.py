from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from sklearn.metrics import roc_auc_score, roc_curve

ROOT_DIR = Path(__file__).resolve().parent.parent
APP_DIR = ROOT_DIR / "app"
RESULTS_DIR = ROOT_DIR / "results" / "final_category_rebuild"
PAPER_DIR = ROOT_DIR / "paper"
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

PUBLISHED_QEP_PARAMS = {
    "QEH": {
        "MW": (7.077e001, 2.830e002, 8.497e001, -1.185e000, 69.5849922),
        "LogP": (9.381e001, 3.077e000, 1.434e000, 6.164e-001, 94.4228257),
        "HBA": (1.176e002, 2.409e000, 1.567e000, 7.155e000, 120.4572352),
        "HBD": (2.334e002, 4.535e-001, -1.480e000, 4.470e000, 228.1589796),
        "RB": (8.470e001, 4.758e000, -2.423e000, 5.437e000, 89.7012502),
        "arR": (3.018e002, 1.101e000, 8.869e-001, -2.281e001, 276.9634213),
    },
    "QEI": {
        "MW": (7.638e001, 2.983e002, 8.364e001, 1.912e000, 78.2919965),
        "LogP": (7.427e001, 4.555e000, -2.193e000, -2.987e000, 71.2829691),
        "HBA": (1.394e002, 1.363e000, 1.283e000, 5.341e-001, 133.9224801),
        "HBD": (6.706e002, -1.163e000, 7.856e-001, 7.951e-001, 331.170104),
        "RB": (6.549e001, 6.219e000, -2.448e000, 5.318e000, 70.5540709),
        "arR": (2.875e002, 3.050e-001, 1.554e000, -8.864e001, 193.0023343),
    },
    "QEF": {
        "MW": (5.103e001, 3.142e002, -5.631e001, 2.342e000, 53.3719946),
        "LogP": (5.073e001, 3.674e000, -1.238e000, 2.067e000, 52.773116),
        "HBA": (7.379e001, 1.841e000, 1.326e000, 5.158e-001, 73.7976536),
        "HBD": (1.647e002, -9.762e-001, -2.027e000, 1.384e000, 144.9887053),
        "RB": (4.091e001, 1.822e000, 2.582e000, 6.235e-001, 41.4385926),
        "arR": (1.344e002, 8.383e-001, 1.347e000, -3.117e001, 102.3024319),
    },
}


def _read_smiles(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return [str(row.get("SMILES", "")).strip() for row in reader if str(row.get("SMILES", "")).strip()]


def _tanimoto_filtered_negatives(positive_smiles: list[str], negative_source_smiles: list[str], threshold: float = 0.3):
    positive_mols = [(smiles, Chem.MolFromSmiles(smiles)) for smiles in positive_smiles]
    positive_mols = [(smiles, mol) for smiles, mol in positive_mols if mol is not None]
    positive_fps = [FP_GEN.GetFingerprint(mol) for _smiles, mol in positive_mols]
    filtered_negatives: list[tuple[str, Chem.Mol]] = []
    for smiles in negative_source_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        fp = FP_GEN.GetFingerprint(mol)
        if max(DataStructs.BulkTanimotoSimilarity(fp, positive_fps)) < threshold:
            filtered_negatives.append((smiles, mol))
    return positive_mols, filtered_negatives


def _published_desirability(x: float, a: float, b: float, c: float, o: float, max_value: float) -> float:
    fitted = o + a * math.exp(-math.exp((-(x - b)) / c) - ((x - b) / c) + 1.0)
    return fitted / max_value


def _published_descriptor_values(mol: Chem.Mol) -> dict[str, float]:
    return {
        "MW": float(Descriptors.MolWt(mol)),
        "LogP": float(Crippen.MolLogP(mol)),
        "HBA": float(rdMolDescriptors.CalcNumHBA(mol)),
        "HBD": float(rdMolDescriptors.CalcNumHBD(mol)),
        "RB": float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "arR": float(rdMolDescriptors.CalcNumAromaticRings(mol)),
    }


def _qex(mol: Chem.Mol, score_name: str) -> float:
    descriptor_values = _published_descriptor_values(mol)
    log_terms: list[float] = []
    for feature in ("MW", "LogP", "HBA", "HBD", "RB", "arR"):
        desirability = _published_desirability(descriptor_values[feature], *PUBLISHED_QEP_PARAMS[score_name][feature])
        if desirability <= 0 or not math.isfinite(desirability):
            return 0.0
        log_terms.append(math.log(desirability))
    return math.exp(sum(log_terms) / len(log_terms))


def _published_pesticide_scores(mol: Chem.Mol) -> dict[str, float]:
    qeh = _qex(mol, "QEH")
    qei = _qex(mol, "QEI")
    qef = _qex(mol, "QEF")
    return {
        "QEH": qeh,
        "QEI": qei,
        "QEF": qef,
        "QEP_max": max(qeh, qei, qef),
        "QEP_avg": (qeh + qei + qef) / 3.0,
    }


def _metrics(pos_scores: list[float], neg_scores: list[float]) -> dict[str, float]:
    labels = [1] * len(pos_scores) + [0] * len(neg_scores)
    scores = pos_scores + neg_scores
    auc = float(roc_auc_score(labels, scores))
    fpr, tpr, thresholds = roc_curve(labels, scores)
    balanced_acc = 0.5 * (tpr + (1.0 - fpr))
    best_idx = int(balanced_acc.argmax())
    return {
        "auc": auc,
        "balanced_accuracy": float(balanced_acc[best_idx]),
        "threshold": float(thresholds[best_idx]),
    }


def main() -> None:
    final_summary_path = RESULTS_DIR / "final_category_rebuild_summary.csv"
    with final_summary_path.open("r", encoding="utf-8-sig", newline="") as handle:
        summary_rows = list(csv.DictReader(handle))
    pesticide_summary = next(row for row in summary_rows if row["category"] == "pesticides")

    positive_smiles = _read_smiles(APP_DIR / "output" / "final_category_rebuild" / "inputs" / "pesticides__positive.csv")
    negative_source_smiles = _read_smiles(APP_DIR / "output" / "final_category_rebuild" / "inputs" / "pesticides__negative_source.csv")
    positive_mols, filtered_negative_mols = _tanimoto_filtered_negatives(positive_smiles, negative_source_smiles, threshold=0.3)

    positive_score_rows = [_published_pesticide_scores(mol) for _smiles, mol in positive_mols]
    negative_score_rows = [_published_pesticide_scores(mol) for _smiles, mol in filtered_negative_mols]

    final_auc = float(pesticide_summary["auc"])
    final_bal_acc = float(pesticide_summary["balanced_accuracy"])

    rows: list[dict[str, object]] = []
    for comparator in ("QEH", "QEI", "QEF", "QEP_max", "QEP_avg"):
        metric = _metrics([row[comparator] for row in positive_score_rows], [row[comparator] for row in negative_score_rows])
        rows.append(
            {
                "category": "pesticides",
                "comparison_regime": "same_pubchem_positive_retained_cross_category_negative_tanimoto_0.3",
                "published_score": comparator,
                "positive_count": len(positive_mols),
                "negative_count": len(filtered_negative_mols),
                "published_auc": round(metric["auc"], 4),
                "published_balanced_accuracy": round(metric["balanced_accuracy"], 4),
                "published_threshold": round(metric["threshold"], 6),
                "final_rebuild_auc": round(final_auc, 4),
                "final_rebuild_balanced_accuracy": round(final_bal_acc, 4),
                "delta_auc_final_minus_published": round(final_auc - metric["auc"], 4),
                "delta_balanced_accuracy_final_minus_published": round(final_bal_acc - metric["balanced_accuracy"], 4),
                "final_beats_published_auc": "yes" if final_auc > metric["auc"] else "no",
                "final_beats_published_balanced_accuracy": "yes" if final_bal_acc > metric["balanced_accuracy"] else "no",
                "implementation_note": "Published QEP coefficients from Avram 2014 supplementary Table S2 reimplemented with RDKit descriptor calculation for same-regime comparison.",
            }
        )

    for output_dir in (RESULTS_DIR, PAPER_DIR):
        output_dir.mkdir(parents=True, exist_ok=True)
        csv_path = output_dir / "pesticide_published_score_comparison.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    meta = {
        "direct_published_score_comparators": {
            "all_categories": ["QED"],
            "pesticides": ["QEH", "QEI", "QEF", "QEP_max", "QEP_avg"],
        },
        "no_other_final_category_has_a_direct_published_continuous_score_confirmed_in_this_project_pass": [
            "animal_drugs",
            "cosmetics",
            "flavoring_agents",
            "food_additives",
            "food_contact_substances",
            "fragrances",
            "human_drugs",
            "solvents",
            "surfactants",
        ],
    }
    with (RESULTS_DIR / "published_score_comparator_inventory.json").open("w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
