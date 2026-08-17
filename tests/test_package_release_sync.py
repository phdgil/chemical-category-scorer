from __future__ import annotations

import json
from pathlib import Path

import chemical_category_scorer as scorer
from app.algorithm_score_engine import list_models


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_MODEL_IDS = [
    "han_endocrine_disruptors",
    "final_flavor_fragrance",
    "final_pesticides",
    "final_surfactants",
]
REMOVED_LIBRARY_FUNCTIONS = {
    "animal_drugs",
    "human_drugs",
    "cosmetics",
    "flavoring_agents",
    "food_additives",
    "food_contact_substances",
    "fragrances",
    "solvents",
}


def test_library_and_desktop_manifest_expose_article_panel_only() -> None:
    release = json.loads((ROOT / "app" / "data" / "app_release_config.json").read_text(encoding="utf-8"))

    assert release["release_version"] == "2.1.0"
    assert release["available_models"] == EXPECTED_MODEL_IDS
    assert scorer.available_models() == EXPECTED_MODEL_IDS
    assert [model["model_id"] for model in list_models(public_only=True)] == EXPECTED_MODEL_IDS


def test_python_library_exports_only_reported_category_helpers() -> None:
    assert scorer.__all__ == [
        "LibraryScore",
        "available_models",
        "details_mol",
        "details_smiles",
        "model_labels",
        "score_mol",
        "score_smiles_value",
        "endocrine_disruptors",
        "flavor_fragrance",
        "pesticides",
        "surfactants",
    ]
    assert all(not hasattr(scorer, name) for name in REMOVED_LIBRARY_FUNCTIONS)


def test_every_released_model_has_calibration_and_definition() -> None:
    calibration = json.loads(
        (ROOT / "app" / "data" / "cross_category_calibration.json").read_text(encoding="utf-8")
    )

    assert calibration["release_version"] == "2.1.0"
    assert list(calibration["models"]) == EXPECTED_MODEL_IDS
    for model_id in EXPECTED_MODEL_IDS:
        assert (ROOT / "app" / "data" / "models" / f"{model_id}.json").is_file()
