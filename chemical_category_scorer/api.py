from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from app.algorithm_score_engine import DEFAULT_MODEL_ID, list_models, refresh_model_registry, score_smiles


@dataclass(frozen=True)
class LibraryScore:
    model_id: str
    label: str
    category: str
    score: float
    threshold: float
    property_score: float
    structure_score: float
    decision: str
    valid: bool
    matched_patterns: tuple[str, ...]


def _ensure_mol(mol: Mol | None) -> Mol:
    if mol is None:
        raise ValueError("A valid RDKit Mol is required.")
    return mol


def _result_to_library_score(result) -> LibraryScore:
    return LibraryScore(
        model_id=result.model_id,
        label=result.model_label,
        category=result.category,
        score=float(result.score),
        threshold=float(result.threshold),
        property_score=float(result.property_score),
        structure_score=float(result.structure_score),
        decision=str(result.decision),
        valid=bool(result.valid),
        matched_patterns=tuple(result.matched_patterns),
    )


def available_models() -> list[str]:
    refresh_model_registry()
    return [model["model_id"] for model in list_models(public_only=True)]


def model_labels() -> dict[str, str]:
    refresh_model_registry()
    return {model["model_id"]: model["label"] for model in list_models(public_only=True)}


def details_smiles(smiles: str, model_id: str = DEFAULT_MODEL_ID) -> LibraryScore:
    refresh_model_registry()
    result = score_smiles(smiles, model_id=model_id)
    return _result_to_library_score(result)


def details_mol(mol: Mol, model_id: str = DEFAULT_MODEL_ID) -> LibraryScore:
    smiles = Chem.MolToSmiles(_ensure_mol(mol), canonical=True)
    return details_smiles(smiles, model_id=model_id)


def score_smiles_value(smiles: str, model_id: str = DEFAULT_MODEL_ID) -> float:
    return details_smiles(smiles, model_id=model_id).score


def score_mol(mol: Mol, model_id: str = DEFAULT_MODEL_ID) -> float:
    return details_mol(mol, model_id=model_id).score


def _category_function(model_id: str) -> Callable[[Mol], float]:
    def scorer(mol: Mol) -> float:
        return score_mol(mol, model_id=model_id)

    scorer.__name__ = model_id.replace("final_", "")
    return scorer


animal_drugs = _category_function("final_animal_drugs")
human_drugs = _category_function("final_human_drugs")
cosmetics = _category_function("final_cosmetics")
endocrine_disruptors = _category_function("han_endocrine_disruptors")
flavoring_agents = _category_function("final_flavoring_agents")
food_additives = _category_function("final_food_additives")
food_contact_substances = _category_function("final_food_contact_substances")
fragrances = _category_function("final_fragrances")
pesticides = _category_function("final_pesticides")
solvents = _category_function("final_solvents")
surfactants = _category_function("final_surfactants")
