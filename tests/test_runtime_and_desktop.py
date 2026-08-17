from __future__ import annotations

import csv
import subprocess
import sys

from app.algorithm_score_engine import (
    ALL_MODELS_SENTINEL,
    MODEL_CONFIGS,
    get_model_role,
    list_models,
    score_csv,
    score_smiles,
    score_smiles_all,
)
from app.desktop_app import _clean_label, format_all_model_results


ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
DDT_SMILES = "Clc1ccc(C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl)cc1"


def test_pesticide_network_augmentation_is_released_and_observable() -> None:
    config = MODEL_CONFIGS["final_pesticides"]
    assert config["model_type"] == "network_augmented_choi"
    assert len(config["network_patterns"]) == 15
    assert config["network_fold_consensus_required"] == 3

    result = score_smiles("Clc1ccc(Cl)cc1", "final_pesticides")
    assert result.valid
    assert result.score >= result.threshold
    assert "network:Clc1ccccc1" in result.matched_patterns


def test_desktop_import_and_direct_script_list_models() -> None:
    import app.desktop_app  # noqa: F401

    completed = subprocess.run(
        [sys.executable, "app/desktop_app.py", "--list-models"],
        check=True,
        capture_output=True,
        text=True,
    )
    listed_ids = [line.split("\t", 1)[0] for line in completed.stdout.splitlines() if line.strip()]
    assert listed_ids == [model["model_id"] for model in list_models(public_only=True)]


def test_public_all_sentinel_uses_ordered_release_models_and_keeps_hidden_explicit() -> None:
    models = list_models(public_only=True)
    model_ids = [model["model_id"] for model in models]

    assert model_ids == [
        "han_endocrine_disruptors",
        "final_flavor_fragrance",
        "final_pesticides",
        "final_surfactants",
    ]
    assert len(model_ids) == 4
    assert "final_endocrine_disruptors" not in model_ids
    assert sum(1 for model in models if model["role"] == "product_use") == 3
    assert [model["model_id"] for model in models if model["role"] == "auxiliary_hazard"] == [
        "han_endocrine_disruptors"
    ]

    all_results = score_smiles_all("CCO")
    assert {result.model_id for result in all_results} == set(model_ids)
    assert "final_endocrine_disruptors" not in {result.model_id for result in all_results}

    hidden_result = score_smiles("CCO", "final_endocrine_disruptors")
    assert hidden_result.model_id == "final_endocrine_disruptors"
    assert get_model_role(hidden_result.model_id) == "hidden"


def test_all_model_formatter_lists_product_scores_and_auxiliary_signal_separately() -> None:
    results = score_smiles_all(ASPIRIN_SMILES)
    formatted = format_all_model_results(results)
    product_models = [model for model in list_models(public_only=True) if model["role"] == "product_use"]
    auxiliary_models = [model for model in list_models(public_only=True) if model["role"] == "auxiliary_hazard"]

    assert "Representative product-use evidence:" in formatted
    assert "Best suggestion" not in formatted
    assert "Raw scores and margins are not calibrated probabilities" in formatted
    assert "cross-category interpretation=" in formatted
    assert "Representative product-use evidence: unresolved; no single category-enriched signal." in formatted
    assert "Product-use category scores:" in formatted
    assert "Auxiliary hazard signal:" in formatted
    assert formatted.index("Product-use category scores:") < formatted.index("Auxiliary hazard signal:")
    assert formatted.count("score=") == 4
    assert formatted.count("threshold=") == 4
    assert formatted.count("margin=") == 4
    assert formatted.count("decision=") == 4
    assert formatted.count("patterns=") == 4

    for model in product_models:
        assert _clean_label(model["label"]) in formatted
    for model in auxiliary_models:
        assert _clean_label(model["label"]) in formatted


def test_all_model_batch_uses_product_representative_and_separate_auxiliary_fields(tmp_path) -> None:
    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.csv"
    with input_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SMILES", "name"])
        writer.writeheader()
        writer.writerow({"SMILES": DDT_SMILES, "name": "ddt"})
        writer.writerow({"SMILES": "not-a-smiles", "name": "invalid"})

    score_csv(input_path, output_path, "SMILES", [ALL_MODELS_SENTINEL])

    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert rows[0]["algorithm_category"] == "pesticides"
    assert rows[0]["algorithm_category"] != "endocrine disruptors"
    assert rows[0]["representative_product_status"] in {"pesticides", "unresolved"}
    assert rows[0]["product_high_specificity_count"].isdigit()
    assert rows[0]["algorithm_margin"]
    assert rows[0]["product_positive_count"] == "1"
    assert rows[0]["auxiliary_hazard_model_id"] == "han_endocrine_disruptors"
    assert rows[0]["auxiliary_hazard_category"] == "endocrine disruptors"
    assert rows[0]["auxiliary_hazard_score"]
    assert rows[0]["auxiliary_hazard_threshold"]
    assert rows[0]["auxiliary_hazard_margin"]
    assert rows[0]["auxiliary_hazard_valid"] == "true"

    assert rows[1]["algorithm_valid"] == "false"
    assert rows[1]["algorithm_category"] == ""
    assert rows[1]["auxiliary_hazard_valid"] == "false"


def test_multi_explicit_batch_falls_back_when_no_product_model_is_selected(tmp_path) -> None:
    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.csv"
    with input_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SMILES", "name"])
        writer.writeheader()
        writer.writerow({"SMILES": DDT_SMILES, "name": "ddt"})

    score_csv(
        input_path,
        output_path,
        "SMILES",
        ["han_endocrine_disruptors", "final_endocrine_disruptors"],
    )

    with output_path.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))

    assert row["algorithm_valid"] == "true"
    assert row["algorithm_category"] == "endocrine disruptors"
    assert row["algorithm_score"]
    assert row["auxiliary_hazard_model_id"] == "han_endocrine_disruptors"
