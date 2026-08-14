from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from functools import lru_cache
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.warning")


@dataclass(frozen=True)
class StructureKeys:
    standardized_smiles: str
    parent_smiles: str
    parent_inchikey: str
    parent_connectivity_key: str


@lru_cache(maxsize=None)
def structure_keys(smiles: str) -> StructureKeys | None:
    """Return deterministic whole-structure and parent-structure identity keys."""
    molecule = Chem.MolFromSmiles((smiles or "").strip())
    if molecule is None:
        return None

    try:
        cleaned = rdMolStandardize.Cleanup(molecule)
        standardized_smiles = Chem.MolToSmiles(cleaned, isomericSmiles=True)
        parent = rdMolStandardize.FragmentParent(cleaned)
        parent = rdMolStandardize.Uncharger().uncharge(parent)
        parent_smiles = Chem.MolToSmiles(parent, isomericSmiles=True)
        parent_inchikey = Chem.MolToInchiKey(parent)
    except Exception:
        return None

    if not standardized_smiles or not parent_smiles or not parent_inchikey:
        return None

    return StructureKeys(
        standardized_smiles=standardized_smiles,
        parent_smiles=parent_smiles,
        parent_inchikey=parent_inchikey,
        parent_connectivity_key=parent_inchikey.split("-", 1)[0],
    )


@dataclass
class ConstructionIndex:
    exact_smiles: set[str]
    parent_connectivity_keys: set[str]
    unresolved_rows: int = 0

    @classmethod
    def empty(cls) -> "ConstructionIndex":
        return cls(set(), set())

    def add(self, keys: StructureKeys | None) -> None:
        if keys is None:
            self.unresolved_rows += 1
            return
        self.exact_smiles.add(keys.standardized_smiles)
        self.parent_connectivity_keys.add(keys.parent_connectivity_key)

    def contains_exact(self, keys: StructureKeys) -> bool:
        return keys.standardized_smiles in self.exact_smiles

    def contains_parent(self, keys: StructureKeys) -> bool:
        return keys.parent_connectivity_key in self.parent_connectivity_keys


def read_delimited(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def load_smiles_index(paths: Iterable[Path], smiles_column: str = "SMILES") -> ConstructionIndex:
    index = ConstructionIndex.empty()
    for path in paths:
        _, rows = read_delimited(path)
        for row in rows:
            index.add(structure_keys(row.get(smiles_column, "")))
    return index


def input_paths(inputs_dir: Path, category: str) -> tuple[Path, Path]:
    positive = inputs_dir / f"{category}__positive.csv"
    negative = inputs_dir / f"{category}__negative_source.csv"
    for path in (positive, negative):
        if not path.is_file():
            raise FileNotFoundError(path)
    return positive, negative


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_external_rows(
    external_rows: list[dict[str, str]],
    *,
    smiles_column: str,
    identifier_column: str,
    target_positive: ConstructionIndex,
    target_negative: ConstructionIndex,
    any_construction: ConstructionIndex,
) -> list[dict[str, str]]:
    audited: list[dict[str, str]] = []
    seen_external_parents: set[str] = set()

    for source_row_number, row in enumerate(external_rows, start=2):
        identifier = (row.get(identifier_column) or "").strip()
        source_smiles = (row.get(smiles_column) or "").strip()
        keys = structure_keys(source_smiles)
        result = dict(row)
        result["source_row_number"] = str(source_row_number)
        result["external_identifier"] = identifier
        result["source_smiles"] = source_smiles

        if keys is None:
            result.update(
                {
                    "standardized_smiles": "",
                    "parent_smiles": "",
                    "parent_inchikey": "",
                    "parent_connectivity_key": "",
                    "overlap_basis": "unresolved",
                    "overlap_scope": "unresolved_structure",
                    "external_status": "unresolved",
                    "duplicate_external_parent": "",
                }
            )
            audited.append(result)
            continue

        exact_positive = target_positive.contains_exact(keys)
        parent_positive = target_positive.contains_parent(keys)
        exact_negative = target_negative.contains_exact(keys)
        parent_negative = target_negative.contains_parent(keys)
        exact_any = any_construction.contains_exact(keys)
        parent_any = any_construction.contains_parent(keys)

        if exact_positive or parent_positive:
            scope = "target_positive"
        elif exact_negative or parent_negative:
            scope = "target_negative_source"
        elif exact_any or parent_any:
            scope = "other_category_construction"
        else:
            scope = "none"

        exact_overlap = exact_positive or exact_negative or exact_any
        parent_overlap = parent_positive or parent_negative or parent_any
        if exact_overlap:
            basis = "exact_standardized_structure"
        elif parent_overlap:
            basis = "parent_connectivity"
        else:
            basis = "none"

        duplicate = keys.parent_connectivity_key in seen_external_parents
        seen_external_parents.add(keys.parent_connectivity_key)
        if scope != "none":
            status = "overlap_excluded"
        elif duplicate:
            status = "duplicate_external_excluded"
        else:
            status = "true_external_candidate"

        result.update(
            {
                "standardized_smiles": keys.standardized_smiles,
                "parent_smiles": keys.parent_smiles,
                "parent_inchikey": keys.parent_inchikey,
                "parent_connectivity_key": keys.parent_connectivity_key,
                "overlap_basis": basis,
                "overlap_scope": scope,
                "external_status": status,
                "duplicate_external_parent": "yes" if duplicate else "no",
            }
        )
        audited.append(result)

    return audited


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fieldnames.append(field)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_audit(
    *,
    external_file: Path,
    inputs_dir: Path,
    target_category: str,
    smiles_column: str,
    identifier_column: str,
    output_dir: Path,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    _, external_rows = read_delimited(external_file)
    if external_rows and smiles_column not in external_rows[0]:
        raise ValueError(f"Missing external SMILES column: {smiles_column}")
    if external_rows and identifier_column not in external_rows[0]:
        raise ValueError(f"Missing external identifier column: {identifier_column}")

    target_positive_path, target_negative_path = input_paths(inputs_dir, target_category)
    target_positive = load_smiles_index([target_positive_path])
    target_negative = load_smiles_index([target_negative_path])
    all_construction_paths = sorted(inputs_dir.glob("*__positive.csv")) + sorted(
        inputs_dir.glob("*__negative_source.csv")
    )
    any_construction = load_smiles_index(all_construction_paths)

    audited = audit_external_rows(
        external_rows,
        smiles_column=smiles_column,
        identifier_column=identifier_column,
        target_positive=target_positive,
        target_negative=target_negative,
        any_construction=any_construction,
    )
    counts = Counter(row["external_status"] for row in audited)
    scopes = Counter(row["overlap_scope"] for row in audited)

    audited_path = output_dir / f"{target_category}_external_overlap_audit.csv"
    filtered_path = output_dir / f"{target_category}_true_external_candidates.csv"
    summary_csv_path = output_dir / f"{target_category}_external_overlap_summary.csv"
    summary_json_path = output_dir / f"{target_category}_external_overlap_summary.json"
    write_csv(audited_path, audited)
    write_csv(
        filtered_path,
        [row for row in audited if row["external_status"] == "true_external_candidate"],
    )

    summary_rows = [
        {"measure": "external_rows", "value": str(len(audited))},
        *({"measure": f"status:{key}", "value": str(value)} for key, value in sorted(counts.items())),
        *({"measure": f"overlap_scope:{key}", "value": str(value)} for key, value in sorted(scopes.items())),
    ]
    write_csv(summary_csv_path, summary_rows)

    summary: dict[str, object] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "external_file": str(external_file),
        "external_file_sha256": sha256_file(external_file),
        "inputs_dir": str(inputs_dir),
        "target_category": target_category,
        "smiles_column": smiles_column,
        "identifier_column": identifier_column,
        "external_rows": len(audited),
        "status_counts": dict(sorted(counts.items())),
        "overlap_scope_counts": dict(sorted(scopes.items())),
        "construction_index": {
            "target_positive_exact": len(target_positive.exact_smiles),
            "target_positive_parent": len(target_positive.parent_connectivity_keys),
            "target_negative_exact": len(target_negative.exact_smiles),
            "target_negative_parent": len(target_negative.parent_connectivity_keys),
            "any_construction_exact": len(any_construction.exact_smiles),
            "any_construction_parent": len(any_construction.parent_connectivity_keys),
            "unresolved_construction_rows": any_construction.unresolved_rows,
        },
        "outputs": {
            "audit_csv": str(audited_path),
            "true_external_candidates_csv": str(filtered_path),
            "summary_csv": str(summary_csv_path),
        },
    }
    summary_json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit external molecules against final model-construction inputs.")
    parser.add_argument("--external-file", type=Path, required=True)
    parser.add_argument("--inputs-dir", type=Path, required=True)
    parser.add_argument("--target-category", required=True)
    parser.add_argument("--smiles-column", default="SMILES")
    parser.add_argument("--identifier-column", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_audit(
        external_file=args.external_file,
        inputs_dir=args.inputs_dir,
        target_category=args.target_category,
        smiles_column=args.smiles_column,
        identifier_column=args.identifier_column,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
