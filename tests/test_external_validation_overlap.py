from __future__ import annotations

import csv
import json
from pathlib import Path

from app.audit_external_validation_overlap import run_audit, structure_keys


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_structure_keys_collapses_salt_to_parent() -> None:
    free_base = structure_keys("CN")
    salt = structure_keys("C[NH3+].[Cl-]")
    assert free_base is not None
    assert salt is not None
    assert free_base.parent_connectivity_key == salt.parent_connectivity_key


def test_run_audit_classifies_overlap_and_true_external(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    outputs = tmp_path / "outputs"
    inputs.mkdir()
    write_csv(inputs / "target__positive.csv", [{"SMILES": "CCO"}])
    write_csv(inputs / "target__negative_source.csv", [{"SMILES": "CC(=O)O"}])
    write_csv(inputs / "other__positive.csv", [{"SMILES": "c1ccccc1"}])
    write_csv(inputs / "other__negative_source.csv", [{"SMILES": "CCN"}])
    external = tmp_path / "external.tsv"
    with external.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ID", "SMILES"], delimiter="\t")
        writer.writeheader()
        writer.writerows(
            [
                {"ID": "positive", "SMILES": "CCO"},
                {"ID": "negative", "SMILES": "CC(=O)O"},
                {"ID": "other", "SMILES": "c1ccccc1"},
                {"ID": "external", "SMILES": "CCCC"},
                {"ID": "duplicate", "SMILES": "CCCC"},
                {"ID": "invalid", "SMILES": "not-a-smiles"},
            ]
        )

    summary = run_audit(
        external_file=external,
        inputs_dir=inputs,
        target_category="target",
        smiles_column="SMILES",
        identifier_column="ID",
        output_dir=outputs,
    )

    assert summary["status_counts"] == {
        "duplicate_external_excluded": 1,
        "overlap_excluded": 3,
        "true_external_candidate": 1,
        "unresolved": 1,
    }
    with (outputs / "target_external_overlap_audit.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        rows = {row["ID"]: row for row in csv.DictReader(handle)}
    assert rows["positive"]["overlap_scope"] == "target_positive"
    assert rows["negative"]["overlap_scope"] == "target_negative_source"
    assert rows["other"]["overlap_scope"] == "other_category_construction"
    assert rows["external"]["external_status"] == "true_external_candidate"
    assert rows["duplicate"]["external_status"] == "duplicate_external_excluded"
    assert rows["invalid"]["external_status"] == "unresolved"
    assert json.loads((outputs / "target_external_overlap_summary.json").read_text())["external_rows"] == 6
