from __future__ import annotations

import csv
import json
from pathlib import Path

from app.resolve_external_candidate_structures import resolve_candidate_file


def test_resolve_candidate_file_uses_cached_structure(tmp_path: Path) -> None:
    input_path = tmp_path / "candidates.csv"
    output_path = tmp_path / "resolved.csv"
    cache_path = tmp_path / "cache.json"
    with input_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["external_id", "external_name", "cas"])
        writer.writeheader()
        writer.writerow({"external_id": "64-17-5", "external_name": "Ethanol", "cas": "64-17-5"})
        writer.writerow({"external_id": "missing", "external_name": "", "cas": ""})
    cache_path.write_text(
        json.dumps(
            {
                "cas:64-17-5": {
                    "status": "resolved",
                    "query": "64-17-5",
                    "cid": 702,
                    "smiles": "CCO",
                    "connectivity_smiles": "CCO",
                    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                    "iupac_name": "ethanol",
                }
            }
        ),
        encoding="utf-8",
    )

    counts = resolve_candidate_file(input_path, output_path, cache_path, checkpoint_every=1)

    assert counts == {"resolved": 1, "missing_identifier": 1}
    with output_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["SMILES"] == "CCO"
    assert rows[0]["resolution_method"] == "cas"
    assert rows[1]["resolution_status"] == "missing_identifier"


def test_resolve_candidate_file_cache_only_materializes_partial_results(tmp_path: Path) -> None:
    input_path = tmp_path / "candidates.csv"
    output_path = tmp_path / "resolved.csv"
    cache_path = tmp_path / "cache.json"
    with input_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["external_id", "external_name", "cas"])
        writer.writeheader()
        writer.writerow({"external_id": "known", "external_name": "Ethanol", "cas": "64-17-5"})
        writer.writerow({"external_id": "uncached", "external_name": "Unknown", "cas": "111-11-1"})
    cache_path.write_text(
        json.dumps(
            {
                "cas:64-17-5": {
                    "status": "resolved",
                    "query": "64-17-5",
                    "smiles": "CCO",
                    "resolution_source": "NCI_CIR",
                }
            }
        ),
        encoding="utf-8",
    )

    counts = resolve_candidate_file(
        input_path,
        output_path,
        cache_path,
        checkpoint_every=1,
        cache_only=True,
    )

    assert counts == {"resolved": 1, "not_attempted": 1}
    with output_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["SMILES"] == "CCO"
    assert rows[1]["resolution_status"] == "not_attempted"
    assert rows[1]["resolution_error"] == "No cached resolution"
