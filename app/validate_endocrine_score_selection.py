from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import QED, rdFingerprintGenerator
from sklearn.metrics import average_precision_score, roc_auc_score

from algorithm_score_engine import refresh_model_registry, score_smiles
from pubchem_category_pipeline import fetch_smiles_for_cids

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
RESULTS_DIR = ROOT_DIR / "results" / "additional_pubchem_checks"
SUMMARY_CSV = RESULTS_DIR / "endocrine_selection_validation.csv"
SUMMARY_JSON = RESULTS_DIR / "endocrine_selection_validation.json"
FINAL_INPUTS_DIR = APP_DIR / "output" / "final_category_rebuild" / "inputs"
HAN_DIR = Path(r"D:/DSWU/2026_기말고사/컴퓨터알고리즘/한세음_20251279_AL")
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

MODELS = {
    "final_endocrine_disruptors": {"label": "GJC uncapped final doctrine"},
    "han_endocrine_disruptors": {"label": "Han Se-eum reconstructed"},
}


def _pick_col(df: pd.DataFrame, cands: list[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for c in cands:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def _read_smiles_csv(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for name in ["SMILES", "Smiles", "smiles", "canonical_smiles", "CanonicalSMILES"]:
            if reader.fieldnames and name in reader.fieldnames:
                return [str(row.get(name, "")).strip() for row in reader if str(row.get(name, "")).strip()]
    return []


def _filter_negatives_after_tanimoto(positive_smiles: list[str], negative_source_smiles: list[str], threshold: float) -> tuple[list[str], list[str]]:
    pos_fps = []
    pos_valid = []
    for smiles in positive_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        pos_fps.append(GEN.GetFingerprint(mol))
        pos_valid.append(smiles)

    neg_fps = []
    neg_valid = []
    for smiles in negative_source_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        neg_fps.append(GEN.GetFingerprint(mol))
        neg_valid.append(smiles)

    if not pos_fps:
        return pos_valid, neg_valid

    max_sim = [0.0] * len(neg_fps)
    for pos_fp in pos_fps:
        sims = DataStructs.BulkTanimotoSimilarity(pos_fp, neg_fps)
        for idx, sim in enumerate(sims):
            if sim > max_sim[idx]:
                max_sim[idx] = float(sim)
    filtered = [smiles for smiles, sim in zip(neg_valid, max_sim) if sim < threshold]
    return pos_valid, filtered


def _load_han_regime() -> tuple[list[str], list[str], float]:
    edc_raw = pd.read_csv(HAN_DIR / "PubChem_EDC.csv", encoding="utf-8-sig")
    smi_col = _pick_col(edc_raw, ["smiles", "canonicalsmiles", "canonical_smiles", "isomericsmiles"])
    if smi_col is None:
        raise ValueError("Could not find EDC smiles column.")
    positive_smiles = [str(v).strip() for v in edc_raw[smi_col].tolist() if str(v).strip()]

    drug_raw = pd.read_csv(HAN_DIR / "PubChem_Drug.csv", encoding="utf-8-sig")
    smi_col = _pick_col(drug_raw, ["smiles", "canonicalsmiles", "canonical_smiles", "isomericsmiles"])
    if smi_col is None:
        raise ValueError("Could not find Drug smiles column.")
    negative_source_smiles = [str(v).strip() for v in drug_raw[smi_col].tolist() if str(v).strip()]
    pos_valid, neg_valid = _filter_negatives_after_tanimoto(positive_smiles, negative_source_smiles, 0.4)
    return pos_valid, neg_valid, 0.4


def _load_final_regime() -> tuple[list[str], list[str], float]:
    positive_smiles = _read_smiles_csv(FINAL_INPUTS_DIR / "endocrine_disruptors__positive.csv")
    negative_source_smiles = _read_smiles_csv(FINAL_INPUTS_DIR / "endocrine_disruptors__negative_source.csv")
    pos_valid, neg_valid = _filter_negatives_after_tanimoto(positive_smiles, negative_source_smiles, 0.3)
    return pos_valid, neg_valid, 0.3


def _evaluate_model(model_id: str, positive_smiles: list[str], negative_smiles: list[str]) -> dict[str, float]:
    pos_results = [score_smiles(smiles, model_id) for smiles in positive_smiles]
    neg_results = [score_smiles(smiles, model_id) for smiles in negative_smiles]
    pos_scores = np.asarray([result.score for result in pos_results if result.valid], dtype=float)
    neg_scores = np.asarray([result.score for result in neg_results if result.valid], dtype=float)
    labels = np.array([1] * len(pos_scores) + [0] * len(neg_scores), dtype=int)
    scores = np.concatenate([pos_scores, neg_scores])
    auc = float(roc_auc_score(labels, scores))
    pr_auc = float(average_precision_score(labels, scores))

    threshold = float(pos_results[0].threshold if pos_results else 0.5)
    pred = scores >= threshold
    tp = int(((labels == 1) & pred).sum())
    fn = int(((labels == 1) & (~pred)).sum())
    tn = int(((labels == 0) & (~pred)).sum())
    fp = int(((labels == 0) & pred).sum())
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)
    ba = float((sens + spec) / 2.0)
    denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc = float(((tp * tn - fp * fn) / math.sqrt(denom)) if denom > 0 else 0.0)
    return {
        "positive_count": int(len(pos_scores)),
        "negative_count": int(len(neg_scores)),
        "auc": auc,
        "pr_auc": pr_auc,
        "balanced_accuracy": ba,
        "mcc_at_model_threshold": mcc,
        "threshold": threshold,
    }


def _evaluate_qed(positive_smiles: list[str], negative_smiles: list[str]) -> dict[str, float]:
    pos_scores = []
    neg_scores = []
    for smiles in positive_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            pos_scores.append(float(QED.qed(mol)))
    for smiles in negative_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            neg_scores.append(float(QED.qed(mol)))
    labels = np.array([1] * len(pos_scores) + [0] * len(neg_scores), dtype=int)
    scores = np.asarray(pos_scores + neg_scores, dtype=float)
    auc = float(roc_auc_score(labels, scores))
    pr_auc = float(average_precision_score(labels, scores))
    threshold = float(np.median(scores))
    pred = scores >= threshold
    tp = int(((labels == 1) & pred).sum())
    fn = int(((labels == 1) & (~pred)).sum())
    tn = int(((labels == 0) & (~pred)).sum())
    fp = int(((labels == 0) & pred).sum())
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)
    ba = float((sens + spec) / 2.0)
    return {
        "positive_count": int(len(pos_scores)),
        "negative_count": int(len(neg_scores)),
        "auc": auc,
        "pr_auc": pr_auc,
        "balanced_accuracy": ba,
        "threshold": threshold,
    }


def main() -> None:
    refresh_model_registry()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    regimes = {
        "han_native_edc_vs_drug": _load_han_regime(),
        "gjc_uncapped_final_regime": _load_final_regime(),
    }
    rows: list[dict[str, Any]] = []
    for regime_name, (positive_smiles, negative_smiles, tanimoto_threshold) in regimes.items():
        qed_metrics = _evaluate_qed(positive_smiles, negative_smiles)
        rows.append({
            "regime": regime_name,
            "model_id": "qed_baseline",
            "label": "RDKit QED",
            "tanimoto_threshold": tanimoto_threshold,
            **qed_metrics,
        })
        for model_id, meta in MODELS.items():
            metrics = _evaluate_model(model_id, positive_smiles, negative_smiles)
            metrics["delta_auc_vs_qed"] = metrics["auc"] - qed_metrics["auc"]
            metrics["delta_balanced_accuracy_vs_qed"] = metrics["balanced_accuracy"] - qed_metrics["balanced_accuracy"]
            rows.append({
                "regime": regime_name,
                "model_id": model_id,
                "label": meta["label"],
                "tanimoto_threshold": tanimoto_threshold,
                **metrics,
            })

    by_regime: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_regime.setdefault(str(row["regime"]), []).append(row)

    decision = {
        "recommended_final_endocrine_model": "final_endocrine_disruptors",
        "rationale": [
            "Han Se-eum's model is the strongest single-regime benchmark on the EDC-vs-drug comparator set.",
            "The uncapped GJC final endocrine scorer remains the aligned cross-regime deployment choice because it is validated under the same broad-category doctrine as the rest of the final app/manuscript set.",
            "Final deployment should favor same-regime comparability across categories rather than a stronger but regime-specific comparator.",
        ],
        "regime_winners": {
            regime: max([row for row in items if row["model_id"] != "qed_baseline"], key=lambda row: (float(row["auc"]), float(row.get("balanced_accuracy", 0.0))))["model_id"]
            for regime, items in by_regime.items()
        },
        "rows": rows,
    }

    SUMMARY_JSON.write_text(json.dumps(decision, indent=2, ensure_ascii=False), encoding="utf-8")
    with SUMMARY_CSV.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "regime",
            "model_id",
            "label",
            "tanimoto_threshold",
            "positive_count",
            "negative_count",
            "auc",
            "pr_auc",
            "balanced_accuracy",
            "mcc_at_model_threshold",
            "threshold",
            "delta_auc_vs_qed",
            "delta_balanced_accuracy_vs_qed",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(SUMMARY_JSON)


if __name__ == "__main__":
    main()
