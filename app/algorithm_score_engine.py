from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold


RDLogger.DisableLog("rdApp.*")

APP_DIR = Path(__file__).resolve().parent
DATA_DIR = APP_DIR / "data"
MODELS_DIR = DATA_DIR / "models"
EVIDENCE_DIR = DATA_DIR / "evidence_panels"
RELEASE_CONFIG_PATH = DATA_DIR / "app_release_config.json"
CROSS_CATEGORY_CALIBRATION_PATH = DATA_DIR / "cross_category_calibration.json"

ALL_MODELS_SENTINEL = "__all__"
DEFAULT_MODEL_ID = "final_pesticides"
PRODUCT_USE_ROLE = "product_use"
AUXILIARY_HAZARD_ROLE = "auxiliary_hazard"
HIDDEN_MODEL_ROLE = "hidden"
EVIDENCE_ONLY_MODEL_IDS = {"full_cosmetics", "full_food_contact_substances"}
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)



@dataclass(frozen=True)
class ScoreBreakdown:
    smiles: str
    model_id: str
    model_label: str
    category: str
    valid: bool
    score: float
    threshold: float
    margin: float
    property_score: float
    structure_score: float
    decision: str
    matched_patterns: tuple[str, ...] = tuple()


@dataclass(frozen=True)
class EvidenceFamilyHit:
    family_id: str
    family_name: str
    similarity: float
    evidence_score: float
    representative_smiles: str
    member_count: int
    matched_motifs: tuple[str, ...]


@dataclass(frozen=True)
class EvidencePanelExplanation:
    panel_id: str
    valid: bool
    query_smiles: str
    family_hits: tuple[EvidenceFamilyHit, ...]
    related_model_support: tuple[ScoreBreakdown, ...]
    notes: tuple[str, ...]


MODEL_CONFIGS: dict[str, dict[str, Any]] = {}
MODEL_RUNTIMES: dict[str, dict[str, Any]] = {}
EVIDENCE_PANELS: dict[str, dict[str, Any]] = {}
PUBLIC_MODEL_IDS: set[str] = set()
PUBLIC_MODEL_ID_ORDER: tuple[str, ...] = tuple()
MODEL_ROLES: dict[str, str] = {}
CROSS_CATEGORY_CALIBRATION: dict[str, dict[str, float]] = {}
_REGISTRY_SIGNATURE: tuple[Any, ...] | None = None


def _slug(text: str) -> str:
    return "_".join((text or "").strip().lower().split())


def _humanize_category(text: str) -> str:
    return (text or "").strip().replace("_", " ")


def _decision_text(category: str, score: float, threshold: float) -> str:
    prefix = "likely" if score >= threshold else "unlikely"
    return f"{prefix} {_humanize_category(category)}"


def _resolve_smiles_column(fieldnames: list[str], preferred: str | None) -> str:
    if preferred and preferred in fieldnames:
        return preferred
    for candidate in ["SMILES", "Smiles", "smiles", "CanonicalSMILES", "canonical_smiles"]:
        if candidate in fieldnames:
            return candidate
    raise KeyError("No SMILES column found.")


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

LEE_PHYSICO_FUNCS = {
    "MolWt": lambda mol: Descriptors.MolWt(mol),
    "MolLogP": lambda mol: Crippen.MolLogP(mol),
    "NumHDonors": lambda mol: rdMolDescriptors.CalcNumHBD(mol),
    "NumHAcceptors": lambda mol: rdMolDescriptors.CalcNumHBA(mol),
    "NumRotatableBonds": lambda mol: rdMolDescriptors.CalcNumRotatableBonds(mol),
    "TPSA": lambda mol: rdMolDescriptors.CalcTPSA(mol),
    "RingCount": lambda mol: rdMolDescriptors.CalcNumRings(mol),
    "NumHeteroatoms": lambda mol: rdMolDescriptors.CalcNumHeteroatoms(mol),
    "NumAromaticRings": lambda mol: rdMolDescriptors.CalcNumAromaticRings(mol),
}


def _compile_scaffold_patterns(smiles_list: list[str]) -> list[Chem.Mol]:
    patterns: list[Chem.Mol] = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        smarts = Chem.MolToSmarts(mol)
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None:
            patterns.append(pattern)
    return patterns


def _compile_smarts_list(smarts_list: list[str]) -> list[Chem.Mol]:
    patterns: list[Chem.Mol] = []
    for smarts in smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None:
            patterns.append(pattern)
    return patterns


def murcko_smiles(mol: Chem.Mol) -> str:
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold, canonical=True) if scaffold and scaffold.GetNumAtoms() > 0 else ""
    except Exception:
        return ""


def _resolve_data_path(path_value: str) -> Path:
    candidate = Path(path_value)
    return candidate if candidate.is_absolute() else DATA_DIR / candidate


def _registry_signature() -> tuple[Any, ...]:
    signature: list[tuple[str, int, int]] = []
    for directory in (MODELS_DIR, EVIDENCE_DIR):
        if directory.exists():
            for path in sorted(directory.glob("*.json")):
                stat = path.stat()
                signature.append((str(path), stat.st_mtime_ns, stat.st_size))
    if RELEASE_CONFIG_PATH.exists():
        stat = RELEASE_CONFIG_PATH.stat()
        signature.append((str(RELEASE_CONFIG_PATH), stat.st_mtime_ns, stat.st_size))
    if CROSS_CATEGORY_CALIBRATION_PATH.exists():
        stat = CROSS_CATEGORY_CALIBRATION_PATH.stat()
        signature.append((str(CROSS_CATEGORY_CALIBRATION_PATH), stat.st_mtime_ns, stat.st_size))
    return tuple(signature)


def _prepare_evidence_panel(config: dict[str, Any]) -> dict[str, Any]:
    runtime = dict(config)
    runtime["families"] = []
    for family in config.get("families", []):
        representative_smiles = str(family.get("representative_smiles", ""))
        rep_mol = Chem.MolFromSmiles(representative_smiles)
        if rep_mol is None:
            continue
        runtime["families"].append(
            {
                **family,
                "representative_fp": FP_GEN.GetFingerprint(rep_mol),
                "motif_patterns": {
                    item["name"]: Chem.MolFromSmarts(item["smarts"])
                    for item in family.get("top_motifs", [])
                    if item.get("smarts") and Chem.MolFromSmarts(item["smarts"]) is not None
                },
            }
        )
    return runtime



def _prepare_runtime_model(config: dict[str, Any]) -> dict[str, Any]:
    runtime: dict[str, Any] = {"config": config}
    model_type = config["model_type"]
    if model_type == "kim_ppv":
        runtime["scaffold_patterns"] = _compile_scaffold_patterns(config["scaffold_smiles"])
        runtime["residue_patterns"] = _compile_smarts_list(config["residue_smarts"])
    elif model_type == "lee_alert_qed":
        runtime["alert_patterns"] = _compile_smarts_list(config["alert_smarts"])
    elif model_type == "choi_auto":
        runtime["selected_patterns"] = {
            name: Chem.MolFromSmarts(smarts)
            for name, smarts in config["selected_patterns"].items()
            if Chem.MolFromSmarts(smarts) is not None
        }
    elif model_type == "han_edc":
        runtime["smarts_patterns"] = {
            name: Chem.MolFromSmarts(smarts)
            for name, smarts in config["smarts_patterns"].items()
            if Chem.MolFromSmarts(smarts) is not None
        }
        reference_smiles_path = _resolve_data_path(str(config["reference_smiles_path"]))
        with reference_smiles_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            reference_smiles = [str(row.get("canonical_smiles", "")).strip() for row in reader]
        reference_mols = [Chem.MolFromSmiles(smiles) for smiles in reference_smiles if smiles]
        reference_mols = [mol for mol in reference_mols if mol is not None]
        runtime["reference_fps_morgan"] = [FP_GEN.GetFingerprint(mol) for mol in reference_mols]
        runtime["reference_fps_maccs"] = [rdMolDescriptors.GetMACCSKeysFingerprint(mol) for mol in reference_mols]
        runtime["reference_fps_rdkit"] = [Chem.RDKFingerprint(mol, fpSize=2048) for mol in reference_mols]
    else:
        raise ValueError(f"Unsupported model type: {model_type}")
    return runtime


def refresh_model_registry(force: bool = False) -> None:
    global MODEL_CONFIGS, MODEL_RUNTIMES, EVIDENCE_PANELS, PUBLIC_MODEL_IDS, PUBLIC_MODEL_ID_ORDER, MODEL_ROLES, CROSS_CATEGORY_CALIBRATION, _REGISTRY_SIGNATURE
    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    signature = _registry_signature()
    if (
        not force
        and _REGISTRY_SIGNATURE == signature
        and MODEL_CONFIGS
        and MODEL_RUNTIMES
    ):
        return

    MODEL_CONFIGS = {}
    MODEL_RUNTIMES = {}
    EVIDENCE_PANELS = {}
    PUBLIC_MODEL_IDS = set()
    PUBLIC_MODEL_ID_ORDER = tuple()
    MODEL_ROLES = {}
    for path in sorted(MODELS_DIR.glob("*.json")):
        with path.open("r", encoding="utf-8") as f:
            config = json.load(f)
        model_id = config["model_id"]
        MODEL_CONFIGS[model_id] = config
        MODEL_RUNTIMES[model_id] = _prepare_runtime_model(config)
    if RELEASE_CONFIG_PATH.exists():
        with RELEASE_CONFIG_PATH.open("r", encoding="utf-8") as f:
            release_config = json.load(f)
        PUBLIC_MODEL_ID_ORDER = tuple(
            model_id for model_id in release_config.get("available_models", []) if model_id in MODEL_CONFIGS
        )
        PUBLIC_MODEL_IDS = set(PUBLIC_MODEL_ID_ORDER)
        configured_roles = release_config.get("model_roles", {})
        MODEL_ROLES = {
            model_id: str(configured_roles.get(model_id, PRODUCT_USE_ROLE))
            for model_id in PUBLIC_MODEL_ID_ORDER
        }
    else:
        PUBLIC_MODEL_IDS = set(MODEL_CONFIGS.keys())
        PUBLIC_MODEL_ID_ORDER = tuple(MODEL_CONFIGS.keys())
        MODEL_ROLES = {model_id: PRODUCT_USE_ROLE for model_id in PUBLIC_MODEL_ID_ORDER}
    CROSS_CATEGORY_CALIBRATION = {}
    if CROSS_CATEGORY_CALIBRATION_PATH.exists():
        with CROSS_CATEGORY_CALIBRATION_PATH.open("r", encoding="utf-8") as f:
            calibration_payload = json.load(f)
        CROSS_CATEGORY_CALIBRATION = {
            str(model_id): dict(values)
            for model_id, values in calibration_payload.get("models", {}).items()
        }
    for path in sorted(EVIDENCE_DIR.glob("*.json")):
        with path.open("r", encoding="utf-8") as f:
            config = json.load(f)
        panel_id = str(config.get("panel_id", path.stem))
        EVIDENCE_PANELS[panel_id] = _prepare_evidence_panel(config)
    _REGISTRY_SIGNATURE = signature


refresh_model_registry()


def list_models(public_only: bool = False) -> list[dict[str, Any]]:
    models = []
    visible_ids = PUBLIC_MODEL_ID_ORDER if public_only and PUBLIC_MODEL_ID_ORDER else tuple(MODEL_CONFIGS.keys())
    for model_id in visible_ids:
        if model_id not in MODEL_CONFIGS:
            continue
        config = MODEL_CONFIGS[model_id]
        models.append(
            {
                "model_id": model_id,
                "label": config["label"],
                "category": config["category"],
                "threshold": config["threshold"],
                "model_type": config["model_type"],
                "role": get_model_role(model_id),
                "source_student": config.get("source_student", ""),
                "description": config.get("description", ""),
            }
        )
    return models


def get_model_role(model_id: str) -> str:
    return MODEL_ROLES.get(model_id, HIDDEN_MODEL_ROLE)


def _resolve_model_ids(model_ids: list[str] | None) -> list[str]:
    if not MODEL_CONFIGS:
        raise RuntimeError("No scoring models found. Run build_scoring_models.py first.")
    if not model_ids:
        if DEFAULT_MODEL_ID in MODEL_CONFIGS:
            return [DEFAULT_MODEL_ID]
        return [next(iter(MODEL_CONFIGS))]
    resolved: list[str] = []
    for model_id in model_ids:
        if model_id == ALL_MODELS_SENTINEL:
            resolved.extend(PUBLIC_MODEL_ID_ORDER or tuple(MODEL_CONFIGS.keys()))
            continue
        if model_id not in MODEL_CONFIGS:
            raise KeyError(f"Unknown model_id: {model_id}")
        resolved.append(model_id)
    seen: set[str] = set()
    unique: list[str] = []
    for model_id in resolved:
        if model_id not in seen:
            unique.append(model_id)
            seen.add(model_id)
    return unique


def _kim_property_score(mol: Chem.Mol, config: dict[str, Any]) -> float:
    try:
        props = {
            "mw": Descriptors.MolWt(mol),
            "xlogp": Crippen.MolLogP(mol),
            "rotbonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        }
    except Exception:
        return 0.0

    hist_models = config["hist_models"]
    num_bins = 30
    total_score = 0.0
    for column in ["mw", "xlogp", "rotbonds", "aromatic_rings"]:
        value = props[column]
        bins = hist_models["bins"][column]
        bin_index = int(np.digitize(value, bins) - 1)
        bin_index = max(0, min(bin_index, num_bins - 1))
        p_pos = hist_models["pos"][column][bin_index]
        p_neg = hist_models["neg"][column][bin_index]
        total_score += p_pos / (p_pos + p_neg) if (p_pos + p_neg) > 1e-10 else 0.0
    return total_score / 4.0


def _kim_structure_score(mol: Chem.Mol, runtime: dict[str, Any], config: dict[str, Any]) -> float:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        scaffold_score = 0.05
    else:
        scaffold_atoms = scaffold.GetNumHeavyAtoms()
        valid_weights: list[float] = []
        for pattern, weight in zip(runtime["scaffold_patterns"], config["scaffold_weights"]):
            try:
                pattern_atoms = pattern.GetNumHeavyAtoms()
                if scaffold.HasSubstructMatch(pattern) and abs(scaffold_atoms - pattern_atoms) <= 5:
                    valid_weights.append(weight)
            except Exception:
                continue
        scaffold_score = max(valid_weights) if valid_weights else 0.10

    sidechains = Chem.ReplaceCore(mol, scaffold)
    target = sidechains if sidechains is not None else mol
    matched_ppv_sum = 0.0
    for pattern, weight in zip(runtime["residue_patterns"], config["residue_weights"]):
        try:
            if target.HasSubstructMatch(pattern):
                matched_ppv_sum += weight
        except Exception:
            continue
    try:
        scaffold_atoms = scaffold.GetNumHeavyAtoms() if scaffold is not None else 0
        sidechain_atoms = max(0, mol.GetNumHeavyAtoms() - scaffold_atoms)
        ratio = sidechain_atoms / (mol.GetNumHeavyAtoms() + 1e-9)
    except Exception:
        ratio = 0.5
    residue_x = matched_ppv_sum * ratio
    residue_score = 1.0 / (1.0 + np.exp(-5.0 * (residue_x - 0.25)))
    w_scaffold = config["w_Scaffold"]
    return w_scaffold * scaffold_score + (1.0 - w_scaffold) * residue_score


def _score_kim(smiles: str, mol: Chem.Mol, runtime: dict[str, Any], config: dict[str, Any]) -> ScoreBreakdown:
    property_score = _kim_property_score(mol, config)
    structure_score = _kim_structure_score(mol, runtime, config)
    score = config["w_Property"] * property_score + config["w_Structure"] * structure_score
    threshold = float(config["threshold"])
    margin = float(score - threshold)
    decision = _decision_text(config["category"], score, threshold)
    return ScoreBreakdown(
        smiles=smiles,
        model_id=config["model_id"],
        model_label=config["label"],
        category=config["category"],
        valid=True,
        score=float(score),
        threshold=threshold,
        margin=margin,
        property_score=float(property_score),
        structure_score=float(structure_score),
        decision=decision,
        matched_patterns=tuple(),
    )


def _lee_property_component(mol: Chem.Mol, config: dict[str, Any]) -> float:
    score = 1.0
    weight_sum = 0.0
    params = config["params"]
    weights = config["weights"]
    for feature in config["physico_features"]:
        weight = float(weights[feature])
        weight_sum += weight
        value = float(LEE_PHYSICO_FUNCS[feature](mol))
        mu = float(params[feature]["mu"])
        sigma = float(params[feature]["sigma"])
        desirability = np.exp(-0.5 * ((value - mu) / sigma) ** 2)
        score *= desirability ** weight
    return float(score ** (1.0 / weight_sum)) if weight_sum > 0 else 0.0


def _score_lee(smiles: str, mol: Chem.Mol, runtime: dict[str, Any], config: dict[str, Any]) -> ScoreBreakdown:
    property_score = _lee_property_component(mol, config)
    alert_count = 0
    for pattern in runtime["alert_patterns"]:
        try:
            if mol.HasSubstructMatch(pattern):
                alert_count += 1
        except Exception:
            continue
    alert_weight = float(config["weights"]["NegativeAlerts"])
    structure_score = float(config["alert_base"] ** (alert_weight * alert_count))
    combined = float(np.clip(property_score * structure_score, 1e-12, 1.0))
    score = float(combined ** float(config["gamma"]))
    threshold = float(config["threshold"])
    margin = float(score - threshold)
    decision = _decision_text(config["category"], score, threshold)
    return ScoreBreakdown(
        smiles=smiles,
        model_id=config["model_id"],
        model_label=config["label"],
        category=config["category"],
        valid=True,
        score=score,
        threshold=threshold,
        margin=margin,
        property_score=float(property_score),
        structure_score=float(structure_score),
        decision=decision,
        matched_patterns=tuple(),
    )


def _choi_property_component(mol: Chem.Mol, config: dict[str, Any]) -> float:
    if not config["selected_props"]:
        return 0.0
    props = {name: CHOI_PROPERTY_FUNCS[name](mol) for name in config["selected_props"]}
    hits = 0
    for feature in config["selected_props"]:
        low = float(config["ranges"][feature][0])
        high = float(config["ranges"][feature][1])
        if low <= float(props[feature]) <= high:
            hits += 1
    return hits / len(config["selected_props"])


def _choi_structure_component(mol: Chem.Mol, runtime: dict[str, Any], config: dict[str, Any]) -> float:
    weights = config["pattern_weights"]
    total_weight = sum(float(value) for value in weights.values()) or 1.0
    matched = 0.0
    for name, pattern in runtime["selected_patterns"].items():
        try:
            if mol.HasSubstructMatch(pattern):
                matched += float(weights[name])
        except Exception:
            continue
    return matched / total_weight


def _score_choi(smiles: str, mol: Chem.Mol, runtime: dict[str, Any], config: dict[str, Any]) -> ScoreBreakdown:
    property_score = _choi_property_component(mol, config)
    matched_pattern_names: list[str] = []
    weights = config["pattern_weights"]
    total_weight = sum(float(value) for value in weights.values()) or 1.0
    matched_weight = 0.0
    for name, pattern in runtime["selected_patterns"].items():
        try:
            if mol.HasSubstructMatch(pattern):
                matched_pattern_names.append(name)
                matched_weight += float(weights[name])
        except Exception:
            continue
    structure_score = matched_weight / total_weight
    best_w = float(config["best_w"])
    score = best_w * property_score + (1.0 - best_w) * structure_score
    threshold = float(config["threshold"])
    margin = float(score - threshold)
    decision = _decision_text(config["category"], score, threshold)
    return ScoreBreakdown(
        smiles=smiles,
        model_id=config["model_id"],
        model_label=config["label"],
        category=config["category"],
        valid=True,
        score=float(score),
        threshold=threshold,
        margin=margin,
        property_score=float(property_score),
        structure_score=float(structure_score),
        decision=decision,
        matched_patterns=tuple(matched_pattern_names),
    )




def _han_property_component(mol: Chem.Mol, config: dict[str, Any]) -> float:
    stats = config["descriptor_stats"]
    values = {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "RotB": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "ArRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }
    score = 1.0
    for name, value in values.items():
        params = stats[name]
        mu = float(params["mu"])
        sig_up = float(params["sig_up"])
        sig_down = float(params["sig_down"])
        sigma = sig_up if value >= mu else sig_down
        z = (value - mu) / max(sigma, 1e-12)
        score *= float(np.exp(-0.5 * z * z) ** (1.0 / len(values)))
    return float(score)


def _han_similarity_component(mol: Chem.Mol, runtime: dict[str, Any]) -> float:
    ref_morgan = runtime["reference_fps_morgan"]
    if not ref_morgan:
        return 0.0
    query_morgan = FP_GEN.GetFingerprint(mol)
    query_maccs = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
    query_rdkit = Chem.RDKFingerprint(mol, fpSize=2048)
    sim_morgan = max(DataStructs.BulkTanimotoSimilarity(query_morgan, ref_morgan))
    sim_maccs = max(DataStructs.BulkTanimotoSimilarity(query_maccs, runtime["reference_fps_maccs"]))
    sim_rdkit = max(DataStructs.BulkTanimotoSimilarity(query_rdkit, runtime["reference_fps_rdkit"]))
    return float((sim_morgan + sim_maccs + sim_rdkit) / 3.0)


def _weighted_geomean_components(scores: dict[str, float], weights: dict[str, float], eps: float = 1e-6) -> float:
    z = 0.0
    for name, weight in weights.items():
        z += float(weight) * math.log(max(float(scores[name]), eps))
    return float(math.exp(z))


def _score_han(smiles: str, mol: Chem.Mol, runtime: dict[str, Any], config: dict[str, Any]) -> ScoreBreakdown:
    property_score = _han_property_component(mol, config)
    scaffold = murcko_smiles(mol)
    scaffold_score = float(config["scaffold_scores"].get(scaffold, config.get("default_scaffold_score", 0.0) if scaffold else 0.0))
    matched_pattern_names: list[str] = []
    matched_weight = 0.0
    total_smarts_weight = sum(float(value) for value in config["smarts_weights"].values()) or 1.0
    for name, pattern in runtime["smarts_patterns"].items():
        try:
            if mol.HasSubstructMatch(pattern):
                matched_pattern_names.append(name)
                matched_weight += float(config["smarts_weights"][name])
        except Exception:
            continue
    smarts_score = float(min(matched_weight / total_smarts_weight, 1.0))
    similarity_score = _han_similarity_component(mol, runtime)
    score = _weighted_geomean_components(
        {
            "property_score": property_score,
            "scaffold_score": scaffold_score,
            "similarity_score": similarity_score,
            "smarts_score": smarts_score,
        },
        config["weights"],
    )
    structure_weights = {
        "scaffold_score": float(config["weights"]["scaffold_score"]),
        "similarity_score": float(config["weights"]["similarity_score"]),
        "smarts_score": float(config["weights"]["smarts_score"]),
    }
    weight_total = sum(structure_weights.values()) or 1.0
    normalized_structure_weights = {name: value / weight_total for name, value in structure_weights.items()}
    structure_score = _weighted_geomean_components(
        {
            "scaffold_score": scaffold_score,
            "similarity_score": similarity_score,
            "smarts_score": smarts_score,
        },
        normalized_structure_weights,
    )
    threshold = float(config["threshold"])
    margin = float(score - threshold)
    decision = _decision_text(config["category"], score, threshold)
    return ScoreBreakdown(
        smiles=smiles,
        model_id=config["model_id"],
        model_label=config["label"],
        category=config["category"],
        valid=True,
        score=float(score),
        threshold=threshold,
        margin=margin,
        property_score=float(property_score),
        structure_score=float(structure_score),
        decision=decision,
        matched_patterns=tuple(matched_pattern_names),
    )


def _invalid_result(smiles: str, config: dict[str, Any]) -> ScoreBreakdown:
    return ScoreBreakdown(
        smiles=smiles,
        model_id=config["model_id"],
        model_label=config["label"],
        category=config["category"],
        valid=False,
        score=-1e9,
        threshold=float(config["threshold"]),
        margin=-1e9,
        property_score=0.0,
        structure_score=0.0,
        decision="invalid smiles",
        matched_patterns=tuple(),
    )


def score_smiles(smiles: str, model_id: str = DEFAULT_MODEL_ID) -> ScoreBreakdown:
    model_id = _resolve_model_ids([model_id])[0]
    config = MODEL_CONFIGS[model_id]
    runtime = MODEL_RUNTIMES[model_id]
    normalized_smiles = (smiles or "").strip()
    if not normalized_smiles:
        return _invalid_result(smiles, config)
    mol = Chem.MolFromSmiles(normalized_smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return _invalid_result(smiles, config)
    model_type = config["model_type"]
    if model_type == "kim_ppv":
        return _score_kim(smiles, mol, runtime, config)
    if model_type == "lee_alert_qed":
        return _score_lee(smiles, mol, runtime, config)
    if model_type == "choi_auto":
        return _score_choi(smiles, mol, runtime, config)
    if model_type == "han_edc":
        return _score_han(smiles, mol, runtime, config)
    raise ValueError(f"Unsupported model type: {model_type}")


def score_smiles_all(smiles: str, model_ids: list[str] | None = None) -> list[ScoreBreakdown]:
    selected_model_ids = _resolve_model_ids(model_ids or [ALL_MODELS_SENTINEL])
    results = [score_smiles(smiles, model_id=model_id) for model_id in selected_model_ids]
    return sorted(results, key=lambda item: (item.valid, item.score, item.margin), reverse=True)


def choose_best_result(results: list[ScoreBreakdown]) -> ScoreBreakdown | None:
    valid_results = [result for result in results if result.valid]
    return valid_results[0] if valid_results else None


def choose_best_product_result(results: list[ScoreBreakdown]) -> ScoreBreakdown | None:
    valid_product_results = [
        result for result in results if result.valid and get_model_role(result.model_id) == PRODUCT_USE_ROLE
    ]
    return valid_product_results[0] if valid_product_results else None


def cross_category_specificity(result: ScoreBreakdown) -> str:
    calibration = CROSS_CATEGORY_CALIBRATION.get(result.model_id)
    if calibration is None or not result.valid:
        return "unavailable"
    if result.score < result.threshold:
        return "below"
    high_specificity_threshold = float(calibration["high_specificity_threshold"])
    return "high_specificity" if result.score >= high_specificity_threshold else "shared"


def choose_representative_product_result(results: list[ScoreBreakdown]) -> ScoreBreakdown | None:
    specific = [
        result
        for result in results
        if result.valid
        and get_model_role(result.model_id) == PRODUCT_USE_ROLE
        and cross_category_specificity(result) == "high_specificity"
    ]
    return specific[0] if len(specific) == 1 else None


def list_evidence_panels() -> list[dict[str, Any]]:
    panels = []
    for panel_id, config in EVIDENCE_PANELS.items():
        panels.append(
            {
                "panel_id": panel_id,
                "entry_count": int(config.get("entry_count", 0)),
                "related_models": list(config.get("related_models", [])),
                "family_count": len(config.get("families", [])),
                "notes": list(config.get("notes", [])),
            }
        )
    return panels


def explain_evidence_panel(smiles: str, panel_id: str) -> EvidencePanelExplanation:
    panel = EVIDENCE_PANELS.get(panel_id)
    if panel is None:
        raise KeyError(f"Unknown evidence panel: {panel_id}")

    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return EvidencePanelExplanation(
            panel_id=panel_id,
            valid=False,
            query_smiles=smiles,
            family_hits=tuple(),
            related_model_support=tuple(),
            notes=tuple(panel.get("notes", [])),
        )

    query_fp = FP_GEN.GetFingerprint(mol)
    family_hits: list[EvidenceFamilyHit] = []
    for family in panel.get("families", []):
        similarity = float(DataStructs.TanimotoSimilarity(query_fp, family["representative_fp"]))
        matched_motifs = tuple(
            name
            for name, pattern in family.get("motif_patterns", {}).items()
            if pattern is not None and mol.HasSubstructMatch(pattern)
        )
        motif_fraction = len(matched_motifs) / max(1, len(family.get("motif_patterns", {})))
        evidence_score = 0.7 * similarity + 0.3 * motif_fraction
        family_hits.append(
            EvidenceFamilyHit(
                family_id=str(family.get("family_id", "")),
                family_name=str(family.get("family_name", "")),
                similarity=similarity,
                evidence_score=float(evidence_score),
                representative_smiles=str(family.get("representative_smiles", "")),
                member_count=int(family.get("member_count", 0)),
                matched_motifs=matched_motifs,
            )
        )
    family_hits.sort(key=lambda item: (item.evidence_score, item.similarity, item.member_count), reverse=True)

    related_support = [
        score_smiles(smiles, model_id)
        for model_id in panel.get("related_models", [])
        if model_id in MODEL_CONFIGS
    ]
    related_support.sort(key=lambda item: (item.score, item.margin), reverse=True)

    return EvidencePanelExplanation(
        panel_id=panel_id,
        valid=True,
        query_smiles=smiles,
        family_hits=tuple(family_hits[:5]),
        related_model_support=tuple(related_support),
        notes=tuple(panel.get("notes", [])),
    )


def _draw_molecule_png(mol: Chem.Mol, highlight_atoms: list[int] | None = None, width: int = 320, height: int = 220) -> bytes:
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=highlight_atoms or [])
    drawer.FinishDrawing()
    return bytes(drawer.GetDrawingText())


def render_molecule_png(smiles: str, width: int = 320, height: int = 220) -> bytes | None:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return None
    return _draw_molecule_png(mol, width=width, height=height)


def render_pattern_match_png(smiles: str, model_id: str, pattern_name: str, width: int = 180, height: int = 140) -> bytes | None:
    if model_id not in MODEL_CONFIGS:
        return None
    config = MODEL_CONFIGS[model_id]
    runtime = MODEL_RUNTIMES[model_id]
    pattern = runtime.get("selected_patterns", {}).get(pattern_name) or runtime.get("smarts_patterns", {}).get(pattern_name)
    if pattern is None:
        return None
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return None
    match = mol.GetSubstructMatch(pattern)
    if not match:
        return None
    return _draw_molecule_png(mol, highlight_atoms=list(match), width=width, height=height)


def score_csv(
    input_csv: str | Path,
    output_csv: str | Path,
    smiles_column: str | None = None,
    model_ids: list[str] | None = None,
) -> Path:
    requested_all_models = not model_ids or ALL_MODELS_SENTINEL in model_ids
    selected_model_ids = _resolve_model_ids(model_ids)
    input_path = Path(input_csv)
    output_path = Path(output_csv)
    with input_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError("Input CSV has no header.")
        resolved_smiles_column = _resolve_smiles_column(fieldnames, smiles_column)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        result_columns = [
            "algorithm_category",
            "algorithm_score",
            "algorithm_threshold",
            "algorithm_margin",
            "algorithm_property_score",
            "algorithm_structure_score",
            "algorithm_decision",
            "algorithm_valid",
            "product_positive_count",
            "product_high_specificity_count",
            "representative_product_status",
            "auxiliary_hazard_model_id",
            "auxiliary_hazard_category",
            "auxiliary_hazard_score",
            "auxiliary_hazard_threshold",
            "auxiliary_hazard_margin",
            "auxiliary_hazard_decision",
            "auxiliary_hazard_valid",
        ]

        with output_path.open("w", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames + result_columns)
            writer.writeheader()
            for row in reader:
                smiles = row.get(resolved_smiles_column, "")
                results = score_smiles_all(smiles, selected_model_ids)
                result = (
                    results[0]
                    if len(selected_model_ids) == 1 and not requested_all_models
                    else choose_best_product_result(results) or choose_best_result(results)
                )
                auxiliary_result = next(
                    (item for item in results if get_model_role(item.model_id) == AUXILIARY_HAZARD_ROLE),
                    None,
                )
                product_positive_count = sum(
                    1
                    for item in results
                    if item.valid and get_model_role(item.model_id) == PRODUCT_USE_ROLE and item.score >= item.threshold
                )
                high_specificity_products = [
                    item
                    for item in results
                    if item.valid
                    and get_model_role(item.model_id) == PRODUCT_USE_ROLE
                    and cross_category_specificity(item) == "high_specificity"
                ]
                representative_result = choose_representative_product_result(results)
                if result and result.valid:
                    row.update(
                        {
                            "algorithm_category": _humanize_category(result.category),
                            "algorithm_score": f"{result.score:.6f}",
                            "algorithm_threshold": f"{result.threshold:.6f}",
                            "algorithm_margin": f"{result.margin:.6f}",
                            "algorithm_property_score": f"{result.property_score:.6f}",
                            "algorithm_structure_score": f"{result.structure_score:.6f}",
                            "algorithm_decision": result.decision,
                            "algorithm_valid": "true",
                        }
                    )
                else:
                    row.update(
                        {
                            "algorithm_category": "",
                            "algorithm_score": "",
                            "algorithm_threshold": "",
                            "algorithm_margin": "",
                            "algorithm_property_score": "",
                            "algorithm_structure_score": "",
                            "algorithm_decision": "invalid smiles",
                            "algorithm_valid": "false",
                        }
                    )
                row["product_positive_count"] = str(product_positive_count) if results else ""
                row["product_high_specificity_count"] = (
                    str(len(high_specificity_products)) if results else ""
                )
                row["representative_product_status"] = (
                    _humanize_category(representative_result.category)
                    if representative_result is not None
                    else "unresolved"
                )
                if auxiliary_result is not None:
                    row.update(
                        {
                            "auxiliary_hazard_model_id": auxiliary_result.model_id,
                            "auxiliary_hazard_category": _humanize_category(auxiliary_result.category),
                            "auxiliary_hazard_score": f"{auxiliary_result.score:.6f}" if auxiliary_result.valid else "",
                            "auxiliary_hazard_threshold": f"{auxiliary_result.threshold:.6f}",
                            "auxiliary_hazard_margin": f"{auxiliary_result.margin:.6f}" if auxiliary_result.valid else "",
                            "auxiliary_hazard_decision": auxiliary_result.decision,
                            "auxiliary_hazard_valid": "true" if auxiliary_result.valid else "false",
                        }
                    )
                else:
                    row.update(
                        {
                            "auxiliary_hazard_model_id": "",
                            "auxiliary_hazard_category": "",
                            "auxiliary_hazard_score": "",
                            "auxiliary_hazard_threshold": "",
                            "auxiliary_hazard_margin": "",
                            "auxiliary_hazard_decision": "",
                            "auxiliary_hazard_valid": "",
                        }
                    )
                writer.writerow(row)
    return output_path
