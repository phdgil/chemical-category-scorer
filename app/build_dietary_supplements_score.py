from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import QED, rdFingerprintGenerator

from build_scoring_models import build_choi_model
from pubchem_category_pipeline import fetch_cids_for_hnid, fetch_smiles_for_cids
from validate_subtyping_reason import read_smiles, write_smiles_csv

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
FINAL_INPUTS_DIR = APP_DIR / "output" / "final_category_rebuild" / "inputs"
RUN_DIR = APP_DIR / "output" / "additional_pubchem_checks" / "dietary_supplements"
INPUTS_DIR = RUN_DIR / "inputs"
MODELS_DIR = RUN_DIR / "models"
RESULTS_DIR = ROOT_DIR / "results" / "additional_pubchem_checks"
SUMMARY_CSV = RESULTS_DIR / "dietary_supplements_score_summary.csv"
SUMMARY_JSON = RESULTS_DIR / "dietary_supplements_score_summary.json"
MODEL_PATH = MODELS_DIR / "dietary_supplements.json"
POSITIVE_CSV = INPUTS_DIR / "dietary_supplements__positive.csv"
NEGATIVE_CSV = INPUTS_DIR / "dietary_supplements__negative_source.csv"

DIETARY_SUPPLEMENTS_HNID = 18246704
MODEL_ID = "dietary_supplements"
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
FINAL_CATEGORIES = [
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
]


def _collect_negative_source(target_positive: list[str]) -> list[str]:
    positive = set(target_positive)
    negatives: list[str] = []
    seen: set[str] = set()
    for category in FINAL_CATEGORIES:
        path = FINAL_INPUTS_DIR / f"{category}__positive.csv"
        for smiles in read_smiles(path):
            if smiles in positive or smiles in seen:
                continue
            seen.add(smiles)
            negatives.append(smiles)
    return negatives


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

    filtered = []
    max_sim = [0.0] * len(neg_fps)
    for pos_fp in pos_fps:
        sims = DataStructs.BulkTanimotoSimilarity(pos_fp, neg_fps)
        for idx, sim in enumerate(sims):
            if sim > max_sim[idx]:
                max_sim[idx] = float(sim)
    for smiles, sim in zip(neg_valid, max_sim):
        if sim < threshold:
            filtered.append(smiles)
    return pos_valid, filtered


def _score_distribution_metrics(pos_scores: list[float], neg_scores: list[float]) -> dict[str, float]:
    import numpy as np
    from sklearn.metrics import roc_auc_score

    if len(pos_scores) < 2 or len(neg_scores) < 2:
        return {"auc": 0.5, "balanced_accuracy": 0.5}

    labels = np.array([1] * len(pos_scores) + [0] * len(neg_scores), dtype=int)
    scores = np.asarray(pos_scores + neg_scores, dtype=float)
    auc = float(roc_auc_score(labels, scores))

    thresholds = np.unique(np.quantile(scores, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]))
    best_ba = 0.0
    for thr in thresholds:
        pred = scores >= float(thr)
        tp = int(((labels == 1) & pred).sum())
        fn = int(((labels == 1) & (~pred)).sum())
        tn = int(((labels == 0) & (~pred)).sum())
        fp = int(((labels == 0) & pred).sum())
        sens = tp / max(tp + fn, 1)
        spec = tn / max(tn + fp, 1)
        ba = (sens + spec) / 2.0
        if ba > best_ba:
            best_ba = ba
    return {"auc": auc, "balanced_accuracy": float(best_ba)}


def _qed_metrics(positive_smiles: list[str], negative_smiles: list[str]) -> dict[str, float]:
    pos_scores: list[float] = []
    neg_scores: list[float] = []
    for smiles in positive_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            pos_scores.append(float(QED.qed(mol)))
    for smiles in negative_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            neg_scores.append(float(QED.qed(mol)))
    return _score_distribution_metrics(pos_scores, neg_scores)


def main() -> None:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    cids = fetch_cids_for_hnid(DIETARY_SUPPLEMENTS_HNID)
    positive_smiles = fetch_smiles_for_cids(cids)
    negative_source_smiles = _collect_negative_source(positive_smiles)

    write_smiles_csv(POSITIVE_CSV, positive_smiles)
    write_smiles_csv(NEGATIVE_CSV, negative_source_smiles)

    build_choi_model(
        positive_csv=POSITIVE_CSV,
        negative_source_csv=NEGATIVE_CSV,
        model_id=MODEL_ID,
        label="Dietary supplements",
        category="dietary_supplements",
        output_path=MODEL_PATH,
        tanimoto_threshold=0.3,
        ks_threshold=0.1,
        ratio_threshold=3.0,
        use_bayesian_optimization=True,
        bayes_trials=18,
        seed=42,
        category_prior=None,
    )

    model = json.loads(MODEL_PATH.read_text(encoding="utf-8"))
    pos_valid, neg_valid = _filter_negatives_after_tanimoto(positive_smiles, negative_source_smiles, float(model.get("tanimoto_threshold", 0.3)))
    qed = _qed_metrics(pos_valid, neg_valid)
    summary = {
        "model_id": model["model_id"],
        "label": model["label"],
        "hnid": DIETARY_SUPPLEMENTS_HNID,
        "positive_count": len(pos_valid),
        "negative_count": len(neg_valid),
        "selected_prop_count": len(model.get("selected_props", [])),
        "pattern_count": len(model.get("pattern_weights", {})),
        "auc": float(model["metrics"].get("auc", 0.0)),
        "balanced_accuracy": float(model["metrics"].get("balanced_accuracy", 0.0)),
        "threshold": float(model.get("threshold", 0.5)),
        "qed_auc": float(qed["auc"]),
        "qed_balanced_accuracy": float(qed["balanced_accuracy"]),
        "delta_auc_vs_qed": float(model["metrics"].get("auc", 0.0) - qed["auc"]),
        "delta_balanced_accuracy_vs_qed": float(model["metrics"].get("balanced_accuracy", 0.0) - qed["balanced_accuracy"]),
        "positive_csv": str(POSITIVE_CSV),
        "negative_source_csv": str(NEGATIVE_CSV),
        "model_path": str(MODEL_PATH),
        "negative_policy": "Union of final category positives with exact target-overlap removed, followed by Tanimoto filtering at 0.3.",
    }

    fieldnames = list(summary.keys())
    with SUMMARY_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(summary)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(SUMMARY_JSON)


if __name__ == "__main__":
    main()
