from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import optuna
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from scipy.integrate import trapezoid
from scipy.stats import gaussian_kde, ks_2samp
from sklearn.metrics import roc_auc_score, roc_curve




RDLogger.DisableLog("rdApp.*")
optuna.logging.set_verbosity(optuna.logging.WARNING)
warnings.filterwarnings("ignore", category=optuna.exceptions.ExperimentalWarning)


APP_DIR = Path(__file__).resolve().parent
DATA_DIR = APP_DIR / "data"
MODELS_DIR = DATA_DIR / "models"

SOURCE_ROOT = Path(r"D:/DSWU/2026_기말고사/컴퓨터알고리즘")
KIM_DIR = SOURCE_ROOT / "김나연" / "github_repo"
LEE_DIR = SOURCE_ROOT / "20251288_이서윤" / "20251288_이서윤"
CHOI_DIR = SOURCE_ROOT / "20251266_최예빈"

GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

CHOI_PROPERTY_FUNCS = {
    "MW": lambda mol: Descriptors.MolWt(mol),
    "logP": lambda mol: Crippen.MolLogP(mol),
    "HBD": lambda mol: rdMolDescriptors.CalcNumHBD(mol),
    "HBA": lambda mol: rdMolDescriptors.CalcNumHBA(mol),
    "TPSA": lambda mol: rdMolDescriptors.CalcTPSA(mol),
    "RotBonds": lambda mol: rdMolDescriptors.CalcNumRotatableBonds(mol),
    "FCsp3": lambda mol: rdMolDescriptors.CalcFractionCSP3(mol),
    "AromaticRings": lambda mol: rdMolDescriptors.CalcNumAromaticRings(mol),
}

CHOI_CANDIDATE_PATTERNS = {
    "halogen": "[Cl,Br,F,I]",
    "chloroaromatic": "[Cl]c",
    "fluoroaromatic": "[F]c",
    "CF3": "C(F)(F)F",
    "triazole": "n1cncn1",
    "pyridine": "c1ccncc1",
    "imidazole": "c1cnc[nH]1",
    "nitro": "[N+](=O)[O-]",
    "carbamate": "OC(=O)N",
    "urea": "NC(=O)N",
    "ester": "C(=O)OC",
    "aromatic_oh": "Oc1ccccc1",
    "aromatic_nh2": "Nc1ccccc1",
    "sulfonate": "S(=O)(=O)",
    "quaternary_n": "[N+]",
    "aldehyde": "[CH1](=O)",
    "long_chain": "CCCCCC",
    "glycol": "OCCO",
    "benzophenone": "c1ccccc1C(=O)c1ccccc1",
    "cinnamate": "OC(=O)/C=C/c1ccccc1",
}
COMPILED_CHOI_CANDIDATES = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in CHOI_CANDIDATE_PATTERNS.items()
    if Chem.MolFromSmarts(smarts) is not None
}
CHOI_RANGE_OPTIONS = [
    (0.05, 0.95),
    (0.10, 0.90),
    (0.15, 0.85),
    (0.20, 0.80),
    (0.25, 0.75),
]
CHOI_DEFAULT_BAYES_TRIALS = 18
CHOI_CATEGORY_PRIORS: dict[str, dict[str, Any]] = {
    "cosmetics": {
        "source": "data/evidence_panels/cosmetics.json plus native Choi cosmetics candidate patterns",
        "motifs": {
            "long_chain": 6.0,
            "ester": 5.0,
            "glycol": 3.0,
            "sulfonate": 3.0,
            "quaternary_n": 2.0,
            "benzophenone": 2.0,
            "cinnamate": 2.0,
        },
    },
    "cosmetic": {
        "alias_for": "cosmetics",
    },
    "food_contact_substances": {
        "source": "data/evidence_panels/food_contact_substances.json",
        "motifs": {
            "long_chain": 6.0,
            "ester": 5.0,
            "glycol": 3.0,
            "sulfonate": 3.0,
        },
    },
    "human_drugs": {
        "source": "data/evidence_panels/human_drugs.json; motifs used only as broad-category priors, not subtype replacement targets",
        "motifs": {
            "imidazole": 5.0,
            "long_chain": 4.0,
            "aromatic_oh": 3.0,
            "urea": 2.0,
            "quaternary_n": 2.0,
            "ester": 2.0,
            "halogen": 2.0,
        },
    },
    "animal_drugs": {
        "source": "data/evidence_panels/animal_drugs.json; motifs used only as broad-category priors, not subtype replacement targets",
        "motifs": {
            "ester": 4.0,
            "long_chain": 4.0,
            "halogen": 3.0,
            "imidazole": 3.0,
            "aromatic_oh": 2.0,
            "urea": 2.0,
        },
    },
    "fragrances": {
        "source": "data/models/choi_fragrance.json and data/models/pubchem_flavoring_agents.json broad Choi pattern evidence",
        "motifs": {
            "ester": 4.0,
            "aldehyde": 3.0,
            "long_chain": 3.0,
            "benzophenone": 1.0,
            "cinnamate": 1.0,
        },
    },
    "fragrance": {
        "alias_for": "fragrances",
    },
    "surfactants": {
        "source": "data/models/choi_surfactant.json broad Choi pattern evidence",
        "motifs": {
            "long_chain": 5.0,
            "glycol": 4.0,
            "sulfonate": 2.0,
            "quaternary_n": 2.0,
        },
    },
    "surfactant": {
        "alias_for": "surfactants",
    },
    "food_additives": {
        "source": "data/models/pubchem_food_additives.json broad Choi pattern evidence",
        "motifs": {
            "aldehyde": 3.0,
            "cinnamate": 1.0,
            "ester": 1.0,
        },
    },
    "flavoring_agents": {
        "source": "data/models/pubchem_flavoring_agents.json broad Choi pattern evidence",
        "motifs": {
            "ester": 4.0,
            "aldehyde": 3.0,
            "cinnamate": 1.0,
        },
    },
    "lipids": {
        "source": "data/models/pubchem_lipids.json broad Choi pattern evidence",
        "motifs": {
            "long_chain": 5.0,
            "glycol": 4.0,
            "ester": 4.0,
            "cinnamate": 1.0,
        },
    },
    "endocrine_disruptors": {
        "source": "data/models/pubchem_endocrine_disruptors.json broad Choi pattern evidence",
        "motifs": {
            "aromatic_oh": 4.0,
            "benzophenone": 2.0,
            "cinnamate": 1.0,
        },
    },
}

LEE_FIXED_WEIGHTS = {
    "MolWt": 0.101,
    "MolLogP": 0.415,
    "NumHDonors": 3.769,
    "NumHAcceptors": 0.101,
    "NumRotatableBonds": 0.101,
    "TPSA": 0.106,
    "RingCount": 2.442,
    "NumHeteroatoms": 0.918,
    "NumAromaticRings": 2.442,
    "NegativeAlerts": 3.436,
}
LEE_PHYSICO_FEATURES = [
    "MolWt",
    "MolLogP",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "TPSA",
    "RingCount",
    "NumHeteroatoms",
    "NumAromaticRings",
]


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    return value


def _ensure_models_dir() -> None:
    MODELS_DIR.mkdir(parents=True, exist_ok=True)


def _resolve_smiles_column(df: pd.DataFrame) -> str:
    for candidate in ["SMILES", "Smiles", "smiles", "CanonicalSMILES", "canonical_smiles"]:
        if candidate in df.columns:
            return candidate
    raise KeyError("No SMILES column found in CSV.")


def _read_smiles(path: str | Path) -> list[str]:
    df = pd.read_csv(path)
    column = _resolve_smiles_column(df)
    values = [str(value).strip() for value in df[column].dropna().tolist()]
    return [value for value in values if value]


def _write_model_config(config: dict[str, Any], output_path: Path) -> Path:
    _ensure_models_dir()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(_json_ready(config), f, indent=2, ensure_ascii=False)
    return output_path


def build_kim_model(output_path: Path | None = None) -> Path:
    with (KIM_DIR / "best_score_model_config.json").open("r", encoding="utf-8") as f:
        source = json.load(f)
    config = {
        "model_id": "kim_pesticide",
        "label": "Pesticide / Kim Nayeon",
        "category": "pesticide",
        "model_type": "kim_ppv",
        "threshold": 0.4340,
        "source_student": "김나연",
        "source_student_id": "20250786",
        "description": "Notebook-derived PPV histogram plus scaffold/residue scorer.",
        "metrics": {
            "roc_auc": 0.9710,
            "pr_auc": 0.9763,
            "accuracy": 0.9169,
            "balanced_accuracy": 0.9169,
            "precision": 0.9175,
            "recall": 0.9160,
            "mcc": 0.8337,
        },
        "w_Property": source["w_Property"],
        "w_Structure": source["w_Structure"],
        "w_Scaffold": source["w_Scaffold"],
        "scaffold_smiles": source["scaffold_smiles"],
        "residue_smarts": source["residue_smarts"],
        "scaffold_weights": source["scaffold_weights"],
        "residue_weights": source["residue_weights"],
        "hist_models": source["hist_models"],
    }
    return _write_model_config(config, output_path or MODELS_DIR / "kim_pesticide.json")


def _lee_get_features(mol: Chem.Mol, alert_patterns: list[Chem.Mol]) -> dict[str, float]:
    return {
        "MolWt": Descriptors.MolWt(mol),
        "MolLogP": Crippen.MolLogP(mol),
        "NumHDonors": rdMolDescriptors.CalcNumHBD(mol),
        "NumHAcceptors": rdMolDescriptors.CalcNumHBA(mol),
        "NumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "RingCount": rdMolDescriptors.CalcNumRings(mol),
        "NumHeteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "NegativeAlerts": float(sum(1 for pat in alert_patterns if mol.HasSubstructMatch(pat))),
    }


def build_lee_model(output_path: Path | None = None) -> Path:
    agro_df = pd.read_csv(LEE_DIR / "agro_data.csv")
    smiles_column = _resolve_smiles_column(agro_df)
    alert_df = pd.read_csv(LEE_DIR / "negative_alert.csv")
    alert_smarts = [str(value) for value in alert_df["smarts"].dropna().tolist()]
    alert_patterns = [Chem.MolFromSmarts(smarts) for smarts in alert_smarts if Chem.MolFromSmarts(smarts) is not None]

    rows = []
    for smiles in agro_df[smiles_column].dropna().tolist():
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is not None:
            rows.append(_lee_get_features(mol, alert_patterns))
    df_pos = pd.DataFrame(rows)

    params = {}
    for feature in LEE_PHYSICO_FEATURES:
        q_low = df_pos[feature].quantile(0.05)
        q_high = df_pos[feature].quantile(0.95)
        filtered = df_pos.loc[df_pos[feature].between(q_low, q_high), feature].astype(float)
        params[feature] = {
            "mu": float(filtered.mean()),
            "sigma": float(max(filtered.std(ddof=1), 1e-6)),
        }

    config = {
        "model_id": "lee_pesticide",
        "label": "Pesticide / Lee Seoyun",
        "category": "pesticide",
        "model_type": "lee_alert_qed",
        "threshold": 0.320,
        "source_student": "이서윤",
        "source_student_id": "20251288",
        "description": "Gaussian property scorer with negative-alert suppression and gamma calibration.",
        "metrics": {
            "auc": 0.953,
            "ks": 0.879,
            "balanced_accuracy": 0.939,
            "best_threshold": 0.320,
            "objective": 0.459,
        },
        "physico_features": LEE_PHYSICO_FEATURES,
        "params": params,
        "weights": LEE_FIXED_WEIGHTS,
        "alert_base": 0.375,
        "gamma": 0.354,
        "alert_smarts": alert_smarts,
    }
    return _write_model_config(config, output_path or MODELS_DIR / "lee_pesticide.json")


def _choi_get_props(mol: Chem.Mol) -> dict[str, float]:
    return {name: func(mol) for name, func in CHOI_PROPERTY_FUNCS.items()}


def _choi_to_valid_fps(smiles_list: list[str]) -> tuple[list[Any], list[str]]:
    fps = []
    valid = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fps.append(GEN.GetFingerprint(mol))
            valid.append(smiles)
    return fps, valid


def _choi_step1_negative(positive_smiles: list[str], negative_source_smiles: list[str], threshold: float) -> tuple[list[str], list[str]]:
    pos_fps, pos_valid = _choi_to_valid_fps(positive_smiles)
    neg_fps, neg_valid = _choi_to_valid_fps(negative_source_smiles)
    negative = [
        smiles
        for smiles, fp in zip(neg_valid, neg_fps)
        if max(DataStructs.BulkTanimotoSimilarity(fp, pos_fps)) < threshold
    ]
    return pos_valid, negative


def _choi_prepare_rows(smiles_list: list[str]) -> list[dict[str, Any]]:
    rows = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        rows.append(
            {
                "smiles": smiles,
                "props": _choi_get_props(mol),
                "matches": {
                    name: mol.HasSubstructMatch(pattern)
                    for name, pattern in COMPILED_CHOI_CANDIDATES.items()
                },
            }
        )
    return rows


def _choi_calc_props(rows: list[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame([row["props"] for row in rows])


def _choi_step2_select_props(
    pos_props: pd.DataFrame,
    neg_props: pd.DataFrame,
    ks_threshold: float,
    lower_q: float,
    upper_q: float,
) -> tuple[list[str], dict[str, tuple[float, float]]]:
    if pos_props.empty or neg_props.empty:
        return [], {}
    selected = [
        feature
        for feature in pos_props.columns
        if ks_2samp(pos_props[feature].dropna(), neg_props[feature].dropna())[0] >= ks_threshold
    ]
    ranges = {
        feature: (
            float(pos_props[feature].quantile(lower_q)),
            float(pos_props[feature].quantile(upper_q)),
        )
        for feature in selected
    }
    return selected, ranges


def _choi_precompute_pattern_rates(pos_rows: list[dict[str, Any]], neg_rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    rates: dict[str, dict[str, float]] = {}
    pos_total = len(pos_rows) or 1
    neg_total = len(neg_rows) or 1
    for name in COMPILED_CHOI_CANDIDATES:
        pos_rate = sum(1 for row in pos_rows if row["matches"][name]) / pos_total * 100.0
        neg_rate = sum(1 for row in neg_rows if row["matches"][name]) / neg_total * 100.0
        rates[name] = {"pos_rate": float(pos_rate), "neg_rate": float(neg_rate)}
    return rates


def _choi_step3_select_patterns(pattern_rates: dict[str, dict[str, float]], ratio_threshold: float) -> dict[str, float]:
    weights = {}
    for name, rate_info in pattern_rates.items():
        pos_rate = rate_info["pos_rate"]
        neg_rate = rate_info["neg_rate"]
        ratio = pos_rate / neg_rate if neg_rate > 0 else (float("inf") if pos_rate > 0 else 0.0)
        if ratio >= ratio_threshold and pos_rate > 0:
            weights[name] = float(pos_rate - neg_rate)
    return weights


def _resolve_choi_category_prior(category_prior: str | None, category: str) -> dict[str, Any] | None:
    if category_prior in (None, "", "none"):
        return None
    prior_key = category if category_prior == "auto" else str(category_prior)
    prior = CHOI_CATEGORY_PRIORS.get(prior_key)
    if prior and "alias_for" in prior:
        prior = CHOI_CATEGORY_PRIORS.get(str(prior["alias_for"]))
    if prior is None:
        raise ValueError(f"Unsupported Choi category prior: {category_prior!r} for category {category!r}")
    motifs = {
        name: float(weight)
        for name, weight in dict(prior.get("motifs", {})).items()
        if name in COMPILED_CHOI_CANDIDATES
    }
    if not motifs:
        raise ValueError(f"Choi category prior {prior_key!r} has no valid candidate motifs.")
    return {
        "prior_id": prior_key,
        "source": prior.get("source", ""),
        "motifs": motifs,
    }


def _choi_apply_category_priors(
    pattern_weights: dict[str, float],
    category_prior: dict[str, Any] | None,
) -> dict[str, float]:
    if category_prior is None:
        return pattern_weights
    merged = dict(pattern_weights)
    for name, prior_weight in category_prior["motifs"].items():
        merged[name] = max(float(merged.get(name, 0.0)), float(prior_weight))
    return merged


def _choi_two_scores(
    rows: list[dict[str, Any]],
    selected_props: list[str],
    ranges: dict[str, tuple[float, float]],
    weights: dict[str, float],
) -> tuple[np.ndarray, np.ndarray]:
    property_scores = []
    structure_scores = []
    total_weight = sum(weights.values()) or 1.0
    for row in rows:
        props = row["props"]
        prop_score = (
            sum(ranges[feature][0] <= props[feature] <= ranges[feature][1] for feature in selected_props) / len(selected_props)
            if selected_props
            else 0.0
        )
        structure_score = sum(
            weights[name]
            for name, matched in row["matches"].items()
            if name in weights and matched
        ) / total_weight
        property_scores.append(prop_score)
        structure_scores.append(structure_score)
    return np.array(property_scores, dtype=float), np.array(structure_scores, dtype=float)


def _safe_kde(scores: np.ndarray) -> gaussian_kde:
    adjusted = np.asarray(scores, dtype=float)
    if adjusted.size < 2:
        adjusted = np.array([0.0, 1e-8], dtype=float)
    elif np.std(adjusted) < 1e-12:
        adjusted = adjusted + np.linspace(-1e-8, 1e-8, adjusted.size)
    return gaussian_kde(adjusted)


def _score_distribution_metrics(pos_scores: np.ndarray, neg_scores: np.ndarray) -> dict[str, float]:
    if len(pos_scores) < 2 or len(neg_scores) < 2:
        return {
            "objective": 10.0,
            "auc": 0.0,
            "ks": 0.0,
            "balanced_accuracy": 0.0,
            "threshold": 0.5,
            "overlap": 1.0,
            "agro_hit_rate": 0.0,
            "negative_reject_rate": 0.0,
        }

    x_grid = np.linspace(-0.05, 1.05, 2000)
    pdf_a = _safe_kde(pos_scores)(x_grid)
    pdf_n = _safe_kde(neg_scores)(x_grid)
    pdf_a = pdf_a / max(float(trapezoid(pdf_a, x_grid)), 1e-12)
    pdf_n = pdf_n / max(float(trapezoid(pdf_n, x_grid)), 1e-12)
    overlap_area = float(trapezoid(np.minimum(pdf_a, pdf_n), x_grid))

    ks_stat, _ = ks_2samp(pos_scores, neg_scores)
    labels = np.r_[np.ones(len(pos_scores)), np.zeros(len(neg_scores))]
    scores = np.r_[pos_scores, neg_scores]
    try:
        auc = float(roc_auc_score(labels, scores))
    except ValueError:
        auc = 0.0

    fpr, tpr, thresholds = roc_curve(labels, scores)
    balanced_acc = 0.5 * (tpr + (1.0 - fpr))
    best_idx = int(np.argmax(balanced_acc))
    threshold_raw = float(thresholds[best_idx])
    threshold = threshold_raw if np.isfinite(threshold_raw) else 0.5
    threshold_penalty = abs(threshold - 0.5) if np.isfinite(threshold_raw) else 1.0
    best_balanced_acc = float(balanced_acc[best_idx])

    objective = float(
        2.0 * overlap_area
        + (1.0 - auc)
        + (1.0 - best_balanced_acc)
        + 0.5 * (1.0 - float(ks_stat))
        + 0.2 * threshold_penalty
    )
    return {
        "objective": objective,
        "auc": auc,
        "ks": float(ks_stat),
        "balanced_accuracy": best_balanced_acc,
        "threshold": float(threshold),
        "overlap": overlap_area,
        "agro_hit_rate": float(np.mean(pos_scores >= threshold)),
        "negative_reject_rate": float(np.mean(neg_scores < threshold)),
    }


def _choi_metrics_for_w(
    pos_rows: list[dict[str, Any]],
    neg_rows: list[dict[str, Any]],
    selected_props: list[str],
    ranges: dict[str, tuple[float, float]],
    weights: dict[str, float],
    w: float,
) -> dict[str, float]:
    pos_p, pos_s = _choi_two_scores(pos_rows, selected_props, ranges, weights)
    neg_p, neg_s = _choi_two_scores(neg_rows, selected_props, ranges, weights)
    pos_scores = w * pos_p + (1.0 - w) * pos_s
    neg_scores = w * neg_p + (1.0 - w) * neg_s
    metrics = _score_distribution_metrics(pos_scores, neg_scores)
    metrics["w"] = float(w)
    return metrics


def _choi_step4_optimize_w(
    pos_rows: list[dict[str, Any]],
    neg_rows: list[dict[str, Any]],
    selected_props: list[str],
    ranges: dict[str, tuple[float, float]],
    weights: dict[str, float],
) -> dict[str, float]:
    best_metrics: dict[str, float] | None = None
    for w in np.arange(0.0, 1.05, 0.05):
        metrics = _choi_metrics_for_w(pos_rows, neg_rows, selected_props, ranges, weights, float(w))
        if best_metrics is None:
            best_metrics = metrics
            continue
        if metrics["objective"] < best_metrics["objective"] - 1e-12:
            best_metrics = metrics
            continue
        if abs(metrics["objective"] - best_metrics["objective"]) <= 1e-12 and metrics["auc"] > best_metrics["auc"]:
            best_metrics = metrics
    return best_metrics or {
        "w": 0.5,
        "objective": 10.0,
        "auc": 0.0,
        "ks": 0.0,
        "balanced_accuracy": 0.0,
        "threshold": 0.5,
        "overlap": 1.0,
        "agro_hit_rate": 0.0,
        "negative_reject_rate": 0.0,
    }


def _choi_bayesian_optimize(
    pos_rows: list[dict[str, Any]],
    neg_rows: list[dict[str, Any]],
    bayes_trials: int,
    seed: int,
    category_prior: dict[str, Any] | None = None,
) -> dict[str, Any]:
    pos_props = _choi_calc_props(pos_rows)
    neg_props = _choi_calc_props(neg_rows)
    pattern_rates = _choi_precompute_pattern_rates(pos_rows, neg_rows)

    def objective(trial: optuna.trial.Trial) -> float:
        range_idx = trial.suggest_categorical("range_idx", list(range(len(CHOI_RANGE_OPTIONS))))
        lower_q, upper_q = CHOI_RANGE_OPTIONS[range_idx]
        ks_threshold = trial.suggest_float("ks_threshold", 0.05, 0.35)
        ratio_threshold = trial.suggest_float("ratio_threshold", 1.5, 6.0)
        w = trial.suggest_float("w", 0.0, 1.0)

        selected_props, ranges = _choi_step2_select_props(pos_props, neg_props, ks_threshold, lower_q, upper_q)
        pattern_weights = _choi_step3_select_patterns(pattern_rates, ratio_threshold)
        pattern_weights = _choi_apply_category_priors(pattern_weights, category_prior)
        if not selected_props and not pattern_weights:
            trial.set_user_attr("payload", {"objective": 10.0})
            return 10.0

        metrics = _choi_metrics_for_w(pos_rows, neg_rows, selected_props, ranges, pattern_weights, w)
        trial.set_user_attr(
            "payload",
            {
                "selected_props": selected_props,
                "ranges": ranges,
                "pattern_weights": pattern_weights,
                "metrics": metrics,
                "lower_q": lower_q,
                "upper_q": upper_q,
                "ks_threshold": ks_threshold,
                "ratio_threshold": ratio_threshold,
            },
        )
        return float(metrics["objective"])

    study = optuna.create_study(
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=seed, multivariate=True),
    )
    study.optimize(objective, n_trials=max(1, bayes_trials), show_progress_bar=False)
    payload = dict(study.best_trial.user_attrs.get("payload", {}))
    payload["optimization_method"] = "lee_tpe_bayesian"
    payload["bayes_trials"] = int(max(1, bayes_trials))
    return payload


def build_choi_model(
    positive_csv: str | Path,
    negative_source_csv: str | Path,
    model_id: str,
    label: str,
    category: str,
    output_path: Path | None = None,
    tanimoto_threshold: float = 0.3,
    ks_threshold: float = 0.1,
    ratio_threshold: float = 3.0,
    use_bayesian_optimization: bool = True,
    bayes_trials: int = CHOI_DEFAULT_BAYES_TRIALS,
    seed: int = 42,
    category_prior: str | None = None,
) -> Path:
    positive_smiles = _read_smiles(positive_csv)
    negative_source_smiles = _read_smiles(negative_source_csv)

    resolved_category_prior = _resolve_choi_category_prior(category_prior, category)

    pos_valid, negative = _choi_step1_negative(positive_smiles, negative_source_smiles, tanimoto_threshold)
    pos_rows = _choi_prepare_rows(pos_valid)
    neg_rows = _choi_prepare_rows(negative)
    if not pos_rows or not neg_rows:
        raise ValueError(f"Insufficient valid molecules to build {model_id}.")

    pos_props = _choi_calc_props(pos_rows)
    neg_props = _choi_calc_props(neg_rows)
    pattern_rates = _choi_precompute_pattern_rates(pos_rows, neg_rows)
    selected_props, ranges = _choi_step2_select_props(
        pos_props,
        neg_props,
        ks_threshold,
        CHOI_RANGE_OPTIONS[0][0],
        CHOI_RANGE_OPTIONS[0][1],
    )
    pattern_weights = _choi_step3_select_patterns(pattern_rates, ratio_threshold)
    pattern_weights = _choi_apply_category_priors(pattern_weights, resolved_category_prior)
    fallback_metrics = _choi_step4_optimize_w(pos_rows, neg_rows, selected_props, ranges, pattern_weights)
    lower_q, upper_q = CHOI_RANGE_OPTIONS[0]
    tuned_ks_threshold = ks_threshold
    tuned_ratio_threshold = ratio_threshold
    optimization_method = "grid_w_objective"
    tuned_trials = 0
    metrics = fallback_metrics

    if use_bayesian_optimization:
        tuned = _choi_bayesian_optimize(
            pos_rows,
            neg_rows,
            bayes_trials=bayes_trials,
            seed=seed,
            category_prior=resolved_category_prior,
        )
        tuned_metrics = tuned.get("metrics", {})
        tuned_objective = float(tuned_metrics.get("objective", 10.0)) if tuned_metrics else 10.0
        fallback_objective = float(fallback_metrics.get("objective", 10.0))
        if tuned_metrics and tuned_objective <= fallback_objective + 1e-12:
            selected_props = tuned.get("selected_props", [])
            ranges = tuned.get("ranges", {})
            pattern_weights = tuned.get("pattern_weights", {})
            metrics = tuned_metrics
            lower_q = float(tuned.get("lower_q", CHOI_RANGE_OPTIONS[0][0]))
            upper_q = float(tuned.get("upper_q", CHOI_RANGE_OPTIONS[0][1]))
            tuned_ks_threshold = float(tuned.get("ks_threshold", ks_threshold))
            tuned_ratio_threshold = float(tuned.get("ratio_threshold", ratio_threshold))
            optimization_method = str(tuned.get("optimization_method", "lee_tpe_bayesian"))
            tuned_trials = int(tuned.get("bayes_trials", bayes_trials))
        else:
            optimization_method = "grid_w_objective_fallback"

    best_threshold = float(metrics.get("threshold", 0.5))
    config = {
        "model_id": model_id,
        "label": label,
        "category": category,
        "model_type": "choi_auto",
        "threshold": best_threshold,
        "source_student": "최예빈",
        "source_student_id": "20251266",
        "description": "Automatically built from PubChem classification CSV exports using the Choi Yebin pipeline with Lee Seoyun-style Bayesian optimization.",
        "metrics": {
            "auc": float(metrics.get("auc", 0.0)),
            "ks": float(metrics.get("ks", 0.0)),
            "balanced_accuracy": float(metrics.get("balanced_accuracy", 0.0)),
            "overlap": float(metrics.get("overlap", 1.0)),
            "objective": float(metrics.get("objective", 10.0)),
            "negative_count": len(neg_rows),
        },
        "selected_props": selected_props,
        "ranges": ranges,
        "pattern_weights": pattern_weights,
        "selected_patterns": {name: CHOI_CANDIDATE_PATTERNS[name] for name in pattern_weights},
        "best_w": float(metrics.get("w", 0.5)),
        "tanimoto_threshold": tanimoto_threshold,
        "ks_threshold": tuned_ks_threshold,
        "ratio_threshold": tuned_ratio_threshold,
        "range_quantiles": [lower_q, upper_q],
        "optimization_method": optimization_method,
        "optimization_trials": tuned_trials,
        "category_prior_used": resolved_category_prior is not None,
        "category_prior": resolved_category_prior,
        "positive_csv": str(Path(positive_csv)),
        "negative_source_csv": str(Path(negative_source_csv)),
    }
    path = output_path or MODELS_DIR / f"{model_id}.json"
    return _write_model_config(config, path)


def build_default_models() -> list[Path]:
    built = []
    built.append(build_kim_model())
    built.append(build_lee_model())
    built.append(
        build_choi_model(
            positive_csv=CHOI_DIR / "PubChem.csv" / "PubChem_Fragrance.csv",
            negative_source_csv=CHOI_DIR / "PubChem.csv" / "PubChem_Pesticides.csv",
            model_id="choi_fragrance",
            label="Fragrance / Choi Yebin",
            category="fragrance",
        )
    )
    built.append(
        build_choi_model(
            positive_csv=CHOI_DIR / "PubChem.csv" / "PubChem_Surfactant.csv",
            negative_source_csv=CHOI_DIR / "PubChem.csv" / "PubChem_Pesticides.csv",
            model_id="choi_surfactant",
            label="Surfactant / Choi Yebin",
            category="surfactant",
        )
    )
    return built


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build scoring model configs for the desktop app.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("defaults", help="Build the default Kim, Lee, and Choi model set.")

    choi_parser = subparsers.add_parser("choi", help="Build a Choi-style scorer from PubChem classification CSV exports.")
    choi_parser.add_argument("--positive", required=True, help="Positive category CSV exported from PubChem classification browser.")
    choi_parser.add_argument("--negative-source", required=True, help="Negative source CSV used before Tanimoto filtering.")
    choi_parser.add_argument("--model-id", required=True, help="Unique model id, e.g. choi_new_category.")
    choi_parser.add_argument("--label", required=True, help="Display label for the app.")
    choi_parser.add_argument("--category", required=True, help="Category name, e.g. surfactant.")
    choi_parser.add_argument("--output", help="Optional output JSON path.")
    choi_parser.add_argument("--bayes-trials", type=int, default=CHOI_DEFAULT_BAYES_TRIALS, help="Number of Lee-style TPE trials for Choi scorer tuning.")
    choi_parser.add_argument("--seed", type=int, default=42, help="Random seed for Choi scorer tuning.")
    choi_parser.add_argument("--no-bayes", action="store_true", help="Disable Bayesian tuning and use fixed thresholds plus grid w search.")
    choi_parser.add_argument("--category-prior", choices=["auto", "cosmetics", "cosmetic", "food_contact_substances", "none"], default="none", help="Optional broad-category motif prior to merge into Choi structure weights.")


    subparsers.add_parser("lee", help="Build the Lee Seoyun pesticide scorer config from source artifacts.")
    subparsers.add_parser("kim", help="Build the Kim Nayeon pesticide scorer config from source artifacts.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "defaults":
        for path in build_default_models():
            print(path)
        return
    if args.command == "kim":
        print(build_kim_model())
        return
    if args.command == "lee":
        print(build_lee_model())
        return
    if args.command == "choi":
        output_path = Path(args.output) if args.output else None
        print(
            build_choi_model(
                positive_csv=args.positive,
                negative_source_csv=args.negative_source,
                model_id=args.model_id,
                label=args.label,
                category=args.category,
                output_path=output_path,
                use_bayesian_optimization=not args.no_bayes,
                bayes_trials=args.bayes_trials,
                seed=args.seed,
                category_prior=args.category_prior,
            )
        )
        return
    raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
