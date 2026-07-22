from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem

from app import structural_pattern_validation as validation
from app import algorithm_score_engine


def _record(smiles: str, category: str = "target", split: str = "train") -> validation.MoleculeRecord:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    scaffold_key = validation.scaffold_group_key(mol)
    return validation.MoleculeRecord(
        category=category,
        smiles=Chem.MolToSmiles(mol, canonical=True),
        smiles_hash=validation.stable_hash(Chem.MolToSmiles(mol, canonical=True)),
        scaffold_key=scaffold_key,
        scaffold_hash=validation.stable_hash(scaffold_key),
        split=split,
        props=validation.molecule_properties(mol),
        mol=mol,
    )


def test_global_scaffold_split_is_deterministic_across_categories_and_row_order() -> None:
    categories = {
        "a": [
            "COC(=O)c1ccccc1",
            "COC(=O)C1CCCCC1",
            "COC(=O)c1ccncc1",
            "CCOC(=O)C",
        ],
        "b": [
            "Clc1ccccc1",
            "ClC1CCCCC1",
            "Clc1ccncc1",
            "CCOC(=O)C",
        ],
    }
    first = validation.build_split_records(categories, seed=23)
    second = validation.build_split_records(
        {category: list(reversed(values)) for category, values in reversed(list(categories.items()))},
        seed=23,
    )

    def assignments(records):
        return sorted(
            (row.category, row.smiles, row.scaffold_key, row.split)
            for split_map in records.values()
            for rows in split_map.values()
            for row in rows
        )

    assert assignments(first) == assignments(second)
    scaffold_splits: dict[str, set[str]] = {}
    molecule_splits: dict[str, set[str]] = {}
    for split_map in first.values():
        for split, rows in split_map.items():
            for row in rows:
                scaffold_splits.setdefault(row.scaffold_key, set()).add(split)
                molecule_splits.setdefault(row.smiles, set()).add(split)
    assert all(len(splits) == 1 for splits in scaffold_splits.values())
    assert all(len(splits) == 1 for splits in molecule_splits.values())
    assert validation.scaffold_group_key(Chem.MolFromSmiles("CCO")).startswith("acyclic_")


def test_brics_candidates_exclude_parent_dummy_atoms_and_compile_as_runtime_smarts() -> None:
    mol = Chem.MolFromSmiles("CCOC(=O)c1ccccc1")
    assert mol is not None
    parent = Chem.MolToSmiles(mol, canonical=True)
    units = validation.brics_fragment_units(mol, min_fragment_atoms=2, max_fragment_atoms=20)
    assert units
    assert parent not in units
    for unit in units:
        unit_mol = Chem.MolFromSmiles(unit)
        assert unit_mol is not None
        assert all(atom.GetAtomicNum() != 0 for atom in unit_mol.GetAtoms())
        assert Chem.MolFromSmarts(validation._unit_to_smarts(unit)) is not None


def test_candidate_discovery_only_sees_training_positive_records() -> None:
    config = replace(
        validation.ValidationConfig(),
        min_positive_count=1,
        min_positive_prevalence=0.0,
        candidate_pool_size=24,
        max_patterns=8,
    )
    train = [_record("CCOC(=O)c1ccccc1")]
    holdout = [_record("N#Cc1ccncc1", split="test")]
    specs = validation.discover_candidate_specs(train, "brics", config)
    assert specs
    assert all("N#C" not in spec.source_unit and "C#N" not in spec.source_unit for spec in specs)

    holdout_specs = validation.discover_candidate_specs(holdout, "brics", config)
    assert {spec.source_unit for spec in specs}.isdisjoint(
        {spec.source_unit for spec in holdout_specs}
    )


def test_shared_pattern_is_rejected_by_cross_category_specificity_gate() -> None:
    config = replace(
        validation.ValidationConfig(),
        min_positive_count=1,
        min_positive_prevalence=0.0,
        min_enrichment=1.0,
        min_specificity_gap=0.10,
        max_category_breadth=9,
    )
    spec = validation.CandidateSpec(
        pattern_id="manual__ester",
        method="hybrid",
        origin="fixed_smarts",
        source_unit="ester",
        smarts="C(=O)O",
        discovery_count=2,
    )
    positives = [_record("CCOC(=O)C"), _record("COC(=O)C")]
    background = {
        "other": [
            _record("CCOC(=O)C", category="other"),
            _record("COC(=O)C", category="other"),
        ]
    }
    assessments, patterns, weights = validation.assess_and_select_patterns(
        [spec], positives, background, config
    )
    assert not assessments[0].selected
    assert "specificity_gap" in assessments[0].rejection_reason
    assert patterns == {}
    assert weights == {}


def test_frozen_test_threshold_is_not_reoptimized(monkeypatch) -> None:
    positive_rows = [_record("CCO", split="test"), _record("CCCO", split="test")]
    negative_rows = [_record("CCCl", category="other", split="test"), _record("CCBr", category="other", split="test")]
    model = validation.FrozenModel(
        category="target",
        method="fixed_smarts",
        selected_props=(),
        ranges={},
        patterns={},
        pattern_weights={},
        best_w=1.0,
        threshold=0.70,
        validation_selection_score=1.0,
        validation_metrics={},
    )

    def fake_scores(rows, _model):
        if rows is positive_rows:
            return np.asarray([0.65, 0.62])
        return np.asarray([0.20, 0.10])

    monkeypatch.setattr(validation, "score_records", fake_scores)
    metrics, _scores = validation.evaluate_frozen_model(
        model,
        positive_rows,
        {regime: negative_rows for regime in validation.EVALUATION_REGIMES},
    )
    assert metrics["all_other"]["threshold"] == 0.70
    assert metrics["all_other"]["sensitivity"] == 0.0
    assert metrics["all_other"]["balanced_accuracy"] == 0.5
    assert metrics["all_other"]["threshold_source_validation"] == 1.0


def test_generated_pattern_model_matches_runtime_choi_score() -> None:
    row = _record("FC(F)(F)c1ccccc1")
    smarts = validation._unit_to_smarts("FC(F)(F)")
    model = validation.FrozenModel(
        category="pesticides",
        method="brics",
        selected_props=("MW",),
        ranges={"MW": (100.0, 250.0)},
        patterns={"brics__cf3": smarts},
        pattern_weights={"brics__cf3": 2.0},
        best_w=0.25,
        threshold=0.5,
        validation_selection_score=1.0,
        validation_metrics={},
    )
    research_score = float(validation.score_records([row], model)[0])
    config = {
        "model_id": "research_runtime_parity",
        "label": "runtime parity",
        "category": "pesticides",
        "model_type": "choi_auto",
        "threshold": model.threshold,
        "selected_props": list(model.selected_props),
        "ranges": {name: list(bounds) for name, bounds in model.ranges.items()},
        "selected_patterns": model.patterns,
        "pattern_weights": model.pattern_weights,
        "best_w": model.best_w,
    }
    runtime = algorithm_score_engine._prepare_runtime_model(config)
    runtime_result = algorithm_score_engine._score_choi(row.smiles, row.mol, runtime, config)
    assert runtime_result.score == research_score


def test_promotion_gate_refuses_one_seed_weak_and_overlapping_candidate() -> None:
    candidate = {
        regime: {"auc": 0.70, "balanced_accuracy": 0.62}
        for regime in validation.EVALUATION_REGIMES
    }
    fixed = {
        regime: {"auc": 0.72, "balanced_accuracy": 0.64}
        for regime in validation.EVALUATION_REGIMES
    }
    prop = {
        regime: {"auc": 0.71, "balanced_accuracy": 0.63}
        for regime in validation.EVALUATION_REGIMES
    }
    seed_results = [
        {
            "outcomes": {
                "pesticides": {
                    "candidate_method": "hybrid",
                    "overall_validation_winner": "fixed_smarts",
                    "candidate_won_validation": False,
                    "pattern_smarts": ["C(=O)O"],
                    "pattern_count": 1,
                    "best_w": 0.5,
                    "test_metrics": candidate,
                    "baseline_metrics": {"fixed_smarts": fixed, "property_only": prop},
                    "specificity_gap": -0.10,
                    "max_off_target_rate": 0.90,
                }
            },
            "bootstrap_rows": [],
        }
    ]
    config = replace(validation.ValidationConfig(), categories=("pesticides",))
    rows = validation.build_promotion_rows(seed_results, config)
    assert rows[0]["promotion_decision"] == "do_not_promote"
    reasons = rows[0]["failure_reasons"]
    assert "insufficient_scaffold_seeds" in reasons
    assert "candidate_not_validation_winner" in reasons
    assert "cross_category_specificity_gap" in reasons
    assert "auc_delta_vs_fixed" in reasons


def test_output_guard_refuses_mixed_or_unmanifested_runs(tmp_path) -> None:
    validation._guard_output_directory(tmp_path, "run-a")
    stray = tmp_path / "stray.csv"
    stray.write_text("old", encoding="utf-8")
    with pytest.raises(FileExistsError, match="no validation manifest"):
        validation._guard_output_directory(tmp_path, "run-a")

    stray.unlink()
    (tmp_path / "manifest.json").write_text(
        '{"run_signature": "run-b"}\n', encoding="utf-8"
    )
    with pytest.raises(FileExistsError, match="different validation run"):
        validation._guard_output_directory(tmp_path, "run-a")
    validation._guard_output_directory(tmp_path, "run-b")


def test_validation_refuses_deployed_model_manuscript_and_input_output_paths(tmp_path) -> None:
    allowed = replace(validation.ValidationConfig(), output_dir=tmp_path / "external-results")
    allowed.validate()

    protected = (
        validation.APP_DIR / "data" / "models",
        validation.ROOT_DIR / "paper",
        validation.APP_DIR,
        validation.DEFAULT_INPUT_DIR,
        validation.DEFAULT_INPUT_DIR.parent,
    )
    for output_dir in protected:
        config = replace(validation.ValidationConfig(), output_dir=output_dir)
        with pytest.raises(ValueError, match="output_dir|protected"):
            config.validate()
