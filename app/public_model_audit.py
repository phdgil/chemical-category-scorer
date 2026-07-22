from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS
from rdkit.Chem.Scaffolds import MurckoScaffold

from app.algorithm_score_engine import (
    AUXILIARY_HAZARD_ROLE,
    MODEL_CONFIGS,
    PRODUCT_USE_ROLE,
    get_model_role,
    list_models,
    score_smiles,
)


RDLogger.DisableLog("rdApp.*")

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
DOCS_DIR = ROOT_DIR / "docs"
FINAL_REBUILD_INPUTS_DIR = APP_DIR / "output" / "final_category_rebuild" / "inputs"

DEFAULT_PROBE_PANEL = DOCS_DIR / "cross_category_probe_panel.csv"
DEFAULT_AUDIT_SCORES = DOCS_DIR / "cross_category_probe_scores.csv"
DEFAULT_AUDIT_SUMMARY = DOCS_DIR / "cross_category_probe_summary.csv"
DEFAULT_PATTERN_OVERLAP = DOCS_DIR / "public_model_pattern_overlap.csv"
DEFAULT_PATTERN_AUDIT = DOCS_DIR / "cross_category_pattern_audit.md"
APP_TEST_SMILES_PANEL = DOCS_DIR / "app_test_smiles_by_category.csv"

PROBE_PANEL_ROWS = [
    {
        "probe_id": "caffeine",
        "molecule_name": "Caffeine",
        "smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "expected_primary_use": "human_drugs",
        "allowed_overlap_notes": "Stimulant and flavor-related food context; product-use models may also rank human-drug-like chemistry.",
    },
    {
        "probe_id": "aspirin",
        "molecule_name": "Aspirin",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "expected_primary_use": "human_drugs",
        "allowed_overlap_notes": "Drug molecule with ester/aromatic motifs that can overlap food, fragrance, or cosmetic fixed patterns.",
    },
    {
        "probe_id": "ddt",
        "molecule_name": "DDT",
        "smiles": "Clc1ccc(C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl)cc1",
        "expected_primary_use": "pesticides",
        "allowed_overlap_notes": "Organochlorine pesticide; auxiliary endocrine signal is scientifically plausible and is not a product category.",
    },
    {
        "probe_id": "bisphenol_a",
        "molecule_name": "Bisphenol A",
        "smiles": "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1",
        "expected_primary_use": "food_contact_substances",
        "allowed_overlap_notes": "Food-contact/plastic-associated compound with expected auxiliary endocrine-disruption signal.",
    },
    {
        "probe_id": "vanillin",
        "molecule_name": "Vanillin",
        "smiles": "COc1cc(C=O)ccc1O",
        "expected_primary_use": "flavoring_agents",
        "allowed_overlap_notes": "Flavor/fragrance aldehyde and phenol chemistry can trigger multiple product-use thresholds.",
    },
    {
        "probe_id": "sds",
        "molecule_name": "Sodium dodecyl sulfate",
        "smiles": "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]",
        "expected_primary_use": "surfactants",
        "allowed_overlap_notes": "Surfactant with long-chain and sulfonate motifs shared by solvent, fragrance, and food-contact scorers.",
    },
    {
        "probe_id": "ethanol",
        "molecule_name": "Ethanol",
        "smiles": "CCO",
        "expected_primary_use": "solvents",
        "allowed_overlap_notes": "Small polar solvent; weak or no positive calls are possible because fixed-pattern models emphasize broader motifs.",
    },
]


@dataclass(frozen=True)
class CandidateRecord:
    method: str
    category: str
    unit: str
    prevalence_count: int
    prevalence_fraction: float
    category_total: int
    breadth: int
    enrichment: float


def _read_release_model_ids() -> list[str]:
    return [str(model["model_id"]) for model in list_models(public_only=True)]


def public_model_ids() -> list[str]:
    return _read_release_model_ids()


def product_model_ids() -> list[str]:
    return [model_id for model_id in public_model_ids() if get_model_role(model_id) == PRODUCT_USE_ROLE]


def write_default_probe_panel(path: str | Path = DEFAULT_PROBE_PANEL) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(output_path, PROBE_PANEL_ROWS, list(PROBE_PANEL_ROWS[0]))
    return output_path


def _read_probe_panel(path: str | Path) -> list[dict[str, str]]:
    panel_path = Path(path)
    if not panel_path.exists():
        if panel_path == DEFAULT_PROBE_PANEL:
            write_default_probe_panel(panel_path)
        else:
            raise FileNotFoundError(panel_path)
    with panel_path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _fmt(value: float) -> str:
    return f"{value:.6f}"


def build_audit_rows(panel_path: str | Path = DEFAULT_PROBE_PANEL) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    probes = _read_probe_panel(panel_path)
    scores: list[dict[str, str]] = []
    summaries: list[dict[str, str]] = []
    model_ids = public_model_ids()

    for probe in probes:
        results = [score_smiles(probe["smiles"], model_id) for model_id in model_ids]
        product_results = [result for result in results if get_model_role(result.model_id) == PRODUCT_USE_ROLE]
        auxiliary_results = [result for result in results if get_model_role(result.model_id) == AUXILIARY_HAZARD_ROLE]
        ranked_product = sorted(product_results, key=lambda item: (item.valid, item.score, item.margin), reverse=True)
        role_rank_by_model: dict[str, int] = {}
        for role in sorted({get_model_role(result.model_id) for result in results}):
            role_results = [result for result in results if get_model_role(result.model_id) == role]
            role_results.sort(key=lambda item: (item.valid, item.score, item.margin), reverse=True)
            role_rank_by_model.update({result.model_id: rank for rank, result in enumerate(role_results, 1)})
        likely_products = [result for result in product_results if result.valid and result.score >= result.threshold]
        auxiliary = auxiliary_results[0] if auxiliary_results else None

        for release_index, result in enumerate(results, 1):
            scores.append(
                {
                    "probe_id": probe["probe_id"],
                    "molecule_name": probe["molecule_name"],
                    "smiles": probe["smiles"],
                    "expected_primary_use": probe["expected_primary_use"],
                    "release_model_order": str(release_index),
                    "score_rank_within_role": str(role_rank_by_model[result.model_id]),
                    "model_id": result.model_id,
                    "model_label": result.model_label,
                    "model_role": get_model_role(result.model_id),
                    "category": result.category,
                    "valid": str(result.valid).lower(),
                    "score": _fmt(result.score),
                    "threshold": _fmt(result.threshold),
                    "margin": _fmt(result.margin),
                    "decision": result.decision,
                    "matched_patterns": "; ".join(result.matched_patterns),
                }
            )

        top_three = ranked_product[:3]
        summaries.append(
            {
                "probe_id": probe["probe_id"],
                "molecule_name": probe["molecule_name"],
                "smiles": probe["smiles"],
                "expected_primary_use": probe["expected_primary_use"],
                "top_product_model_id": top_three[0].model_id if top_three else "",
                "top_product_category": top_three[0].category if top_three else "",
                "top_product_score": _fmt(top_three[0].score) if top_three else "",
                "top_product_margin": _fmt(top_three[0].margin) if top_three else "",
                "top_three_product_categories": "; ".join(f"{item.category}:{item.score:.6f}" for item in top_three),
                "likely_product_count": str(len(likely_products)),
                "likely_product_categories": "; ".join(result.category for result in sorted(likely_products, key=lambda item: item.category)),
                "auxiliary_endocrine_model_id": auxiliary.model_id if auxiliary else "",
                "auxiliary_endocrine_score": _fmt(auxiliary.score) if auxiliary else "",
                "auxiliary_endocrine_threshold": _fmt(auxiliary.threshold) if auxiliary else "",
                "auxiliary_endocrine_margin": _fmt(auxiliary.margin) if auxiliary else "",
                "auxiliary_endocrine_decision": auxiliary.decision if auxiliary else "",
                "allowed_overlap_notes": probe["allowed_overlap_notes"],
            }
        )
    return scores, summaries


def run_audit(
    panel_path: str | Path = DEFAULT_PROBE_PANEL,
    scores_path: str | Path = DEFAULT_AUDIT_SCORES,
    summary_path: str | Path = DEFAULT_AUDIT_SUMMARY,
) -> tuple[Path, Path]:
    scores, summaries = build_audit_rows(panel_path)
    score_fields = [
        "probe_id",
        "molecule_name",
        "smiles",
        "expected_primary_use",
        "release_model_order",
        "score_rank_within_role",
        "model_id",
        "model_label",
        "model_role",
        "category",
        "valid",
        "score",
        "threshold",
        "margin",
        "decision",
        "matched_patterns",
    ]
    summary_fields = [
        "probe_id",
        "molecule_name",
        "smiles",
        "expected_primary_use",
        "top_product_model_id",
        "top_product_category",
        "top_product_score",
        "top_product_margin",
        "top_three_product_categories",
        "likely_product_count",
        "likely_product_categories",
        "auxiliary_endocrine_model_id",
        "auxiliary_endocrine_score",
        "auxiliary_endocrine_threshold",
        "auxiliary_endocrine_margin",
        "auxiliary_endocrine_decision",
        "allowed_overlap_notes",
    ]
    score_output = Path(scores_path)
    summary_output = Path(summary_path)
    _write_csv(score_output, scores, score_fields)
    _write_csv(summary_output, summaries, summary_fields)
    return score_output, summary_output


def _public_pattern_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for model_id in public_model_ids():
        config = MODEL_CONFIGS[model_id]
        if "selected_patterns" in config:
            pattern_source = "selected_patterns"
            patterns = config["selected_patterns"]
        elif "smarts_patterns" in config:
            pattern_source = "smarts_patterns"
            patterns = config["smarts_patterns"]
        else:
            continue
        for name, smarts in sorted(patterns.items(), key=lambda item: item[0].lower()):
            rows.append(
                {
                    "model_id": model_id,
                    "model_label": str(config["label"]),
                    "model_role": get_model_role(model_id),
                    "category": str(config["category"]),
                    "pattern_source": pattern_source,
                    "pattern_name": str(name),
                    "pattern_smarts": str(smarts),
                    "canonical_smarts": canonical_smarts(str(smarts)),
                }
            )
    return rows


def canonical_smarts(smarts: str) -> str:
    pattern = Chem.MolFromSmarts(smarts)
    return Chem.MolToSmarts(pattern) if pattern is not None else ""


def build_pattern_overlap_rows() -> list[dict[str, str]]:
    pattern_rows = _public_pattern_rows()
    by_canonical: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in pattern_rows:
        by_canonical[row["canonical_smarts"] or row["pattern_smarts"]].append(row)

    output: list[dict[str, str]] = []
    for key in sorted(by_canonical):
        group = sorted(by_canonical[key], key=lambda row: (row["model_role"], row["model_id"], row["pattern_name"]))
        categories = sorted({row["category"] for row in group})
        model_ids = sorted({row["model_id"] for row in group})
        output.append(
            {
                "canonical_smarts": key,
                "pattern_smarts_examples": "; ".join(sorted({row["pattern_smarts"] for row in group})),
                "pattern_names": "; ".join(sorted({row["pattern_name"] for row in group})),
                "model_count": str(len(model_ids)),
                "category_count": str(len(categories)),
                "model_ids": "; ".join(model_ids),
                "categories": "; ".join(categories),
                "roles": "; ".join(sorted({row["model_role"] for row in group})),
                "shared_across_models": str(len(model_ids) > 1).lower(),
            }
        )
    return output


def write_pattern_overlap(
    csv_path: str | Path = DEFAULT_PATTERN_OVERLAP,
    markdown_path: str | Path = DEFAULT_PATTERN_AUDIT,
) -> tuple[Path, Path]:
    rows = build_pattern_overlap_rows()
    fields = [
        "canonical_smarts",
        "pattern_smarts_examples",
        "pattern_names",
        "model_count",
        "category_count",
        "model_ids",
        "categories",
        "roles",
        "shared_across_models",
    ]
    csv_output = Path(csv_path)
    markdown_output = Path(markdown_path)
    _write_csv(csv_output, rows, fields)
    shared = [row for row in rows if row["shared_across_models"] == "true"]
    markdown_output.parent.mkdir(parents=True, exist_ok=True)
    markdown_output.write_text(_pattern_audit_markdown(rows, shared), encoding="utf-8")
    return csv_output, markdown_output


def _pattern_audit_markdown(rows: list[dict[str, str]], shared: list[dict[str, str]]) -> str:
    quality = app_test_smiles_quality(APP_TEST_SMILES_PANEL)
    shared_lines = [
        f"- `{row['canonical_smarts']}` appears in {row['model_count']} models: {row['model_ids']}."
        for row in shared[:20]
    ]
    shared_text = "\n".join(shared_lines) if shared_lines else "- No exact canonical SMARTS reuse was found."
    return (
        "# Cross-category pattern audit\n\n"
        "This report is derived only from the public model JSON configs listed in "
        "`app/data/app_release_config.json`; it does not mutate thresholds or model files.\n\n"
        "The 65-row `docs/app_test_smiles_by_category.csv` file is a structural pattern unit-test panel. "
        "It verifies that configured SMARTS can match fixed probe molecules, but it is not a labeled "
        "classification benchmark and should not be used for performance claims.\n\n"
        f"Public pattern units inspected: {len(rows)}. Exact canonical SMARTS reused across more than one model: {len(shared)}.\n\n"
        "## Pattern-unit panel data quality\n\n"
        f"- Rows: {quality['row_count']}\n"
        f"- Columns: {quality['column_count']}\n"
        f"- Missing cells: {quality['missing_cell_count']}\n"
        f"- Exact duplicate rows: {quality['exact_duplicate_row_count']}\n"
        f"- Unique test SMILES: {quality['unique_test_smiles_count']}\n"
        f"- Invalid test SMILES: {quality['invalid_test_smiles_count']}\n"
        f"- Invalid pattern SMARTS: {quality['invalid_pattern_smarts_count']}\n\n"
        "Repeated probe SMILES are expected in this file because one molecule can be used to test several "
        "pattern units across different public models. These repeats are therefore not classification-label "
        "duplicates.\n\n"
        "## Shared fixed SMARTS\n\n"
        f"{shared_text}\n\n"
        "## Interpretation\n\n"
        "Shared fixed SMARTS can help explain why simple molecules can trigger multiple independent category thresholds. "
        "The seven-probe audit is exploratory cross-category behavior checking, not held-out validation. "
        "A 30-per-category exploratory pilot is only justified when a concrete sampled input panel is recorded; "
        "this repository slice does not add such evidence and makes no accuracy or BRICS-improvement claim.\n"
    )


def app_test_smiles_quality(path: str | Path = APP_TEST_SMILES_PANEL) -> dict[str, int]:
    panel_path = Path(path)
    with panel_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    row_tuples = [tuple((field, row.get(field, "")) for field in fieldnames) for row in rows]
    duplicate_count = len(row_tuples) - len(set(row_tuples))
    missing_count = sum(1 for row in rows for field in fieldnames if str(row.get(field, "")).strip() == "")
    unique_smiles = {
        str(row.get("test_smiles", "")).strip()
        for row in rows
        if str(row.get("test_smiles", "")).strip()
    }
    invalid_smiles = sum(1 for smiles in unique_smiles if Chem.MolFromSmiles(smiles) is None)
    unique_smarts = {
        str(row.get("pattern_smarts", "")).strip()
        for row in rows
        if str(row.get("pattern_smarts", "")).strip()
    }
    invalid_smarts = sum(1 for smarts in unique_smarts if Chem.MolFromSmarts(smarts) is None)
    return {
        "row_count": len(rows),
        "column_count": len(fieldnames),
        "missing_cell_count": missing_count,
        "exact_duplicate_row_count": duplicate_count,
        "unique_test_smiles_count": len(unique_smiles),
        "invalid_test_smiles_count": invalid_smiles,
        "invalid_pattern_smarts_count": invalid_smarts,
    }


def _resolve_smiles_column(fieldnames: Iterable[str], preferred: str | None) -> str:
    names = list(fieldnames)
    if preferred and preferred in names:
        return preferred
    for candidate in ("SMILES", "Smiles", "smiles", "CanonicalSMILES", "canonical_smiles"):
        if candidate in names:
            return candidate
    raise KeyError("No SMILES column found.")


def canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None or mol.GetNumAtoms() == 0:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def read_category_smiles_csv(path: str | Path, category_column: str = "category", smiles_column: str | None = None) -> tuple[dict[str, list[str]], int]:
    categories: dict[str, list[str]] = defaultdict(list)
    invalid_count = 0
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("Input CSV has no header.")
        resolved_smiles = _resolve_smiles_column(reader.fieldnames, smiles_column)
        if category_column not in reader.fieldnames:
            raise KeyError(f"No category column found: {category_column}")
        for row in reader:
            category = str(row.get(category_column, "")).strip()
            smiles = canonical_smiles(str(row.get(resolved_smiles, "")))
            if not category or not smiles:
                invalid_count += 1
                continue
            categories[category].append(smiles)
    return {category: sorted(set(values)) for category, values in categories.items()}, invalid_count


def read_final_rebuild_inputs(limit_per_category: int | None = None) -> tuple[dict[str, list[str]], int]:
    categories: dict[str, list[str]] = {}
    invalid_count = 0
    inputs = [
        (
            str(MODEL_CONFIGS[model_id]["category"]),
            FINAL_REBUILD_INPUTS_DIR / f"{MODEL_CONFIGS[model_id]['category']}__positive.csv",
        )
        for model_id in product_model_ids()
    ]
    missing_paths = [path for _category, path in inputs if not path.exists()]
    if missing_paths:
        missing = ", ".join(str(path) for path in missing_paths)
        raise FileNotFoundError(f"Missing final rebuild input files: {missing}")

    for category, path in inputs:
        values, invalid = _read_smiles_only_csv(path, limit_per_category)
        if not values:
            raise ValueError(f"No valid SMILES found in final rebuild input: {path}")
        categories[category] = values
        invalid_count += invalid
    return categories, invalid_count


def _read_smiles_only_csv(path: Path, limit: int | None = None) -> tuple[list[str], int]:
    values: list[str] = []
    invalid = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return values, invalid
        smiles_col = _resolve_smiles_column(reader.fieldnames, None)
        for row in reader:
            smiles = canonical_smiles(str(row.get(smiles_col, "")))
            if not smiles:
                invalid += 1
                continue
            values.append(smiles)
    unique = sorted(set(values))
    if limit is not None and len(unique) > limit:
        unique = _hash_sample(unique, limit)
    return unique, invalid


def _hash_sample(values: list[str], limit: int) -> list[str]:
    if limit < 0:
        raise ValueError("limit must be non-negative.")
    ranked = sorted((hashlib.sha256(value.encode("utf-8")).hexdigest(), value) for value in values)
    return sorted(value for _digest, value in ranked[:limit])


def fixed_smarts_units(mol: Chem.Mol, pattern_rows: list[dict[str, str]] | None = None) -> set[str]:
    rows = pattern_rows if pattern_rows is not None else _public_pattern_rows()
    units: set[str] = set()
    for row in rows:
        pattern = Chem.MolFromSmarts(row["pattern_smarts"])
        if pattern is not None and mol.HasSubstructMatch(pattern):
            units.add(row["canonical_smarts"] or row["pattern_smarts"])
    return units


def murcko_units(mol: Chem.Mol, min_fragment_atoms: int = 3) -> set[str]:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumHeavyAtoms() < min_fragment_atoms:
        return set()
    value = Chem.MolToSmiles(scaffold, canonical=True)
    return {value} if value else set()


def brics_units(mol: Chem.Mol, min_fragment_atoms: int = 3) -> set[str]:
    units: set[str] = set()
    for fragment in BRICS.BRICSDecompose(mol):
        canonical = _canonical_fragment_smiles(fragment, min_fragment_atoms)
        if canonical:
            units.add(canonical)
    return units


def _canonical_fragment_smiles(fragment: str, min_fragment_atoms: int = 3) -> str:
    mol = Chem.MolFromSmiles(fragment)
    if mol is None:
        return ""
    editable = Chem.RWMol(mol)
    for atom_index in sorted((atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0), reverse=True):
        editable.RemoveAtom(atom_index)
    clean = editable.GetMol()
    try:
        Chem.SanitizeMol(clean)
    except Exception:
        return ""
    if clean.GetNumHeavyAtoms() < min_fragment_atoms:
        return ""
    return Chem.MolToSmiles(clean, canonical=True)


def extract_candidate_units(smiles: str, method: str, min_fragment_atoms: int = 3) -> set[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return set()
    if method == "fixed_smarts":
        return fixed_smarts_units(mol)
    if method == "murcko":
        return murcko_units(mol, min_fragment_atoms)
    if method == "brics":
        return brics_units(mol, min_fragment_atoms)
    if method == "hybrid":
        return fixed_smarts_units(mol) | murcko_units(mol, min_fragment_atoms) | brics_units(mol, min_fragment_atoms)
    raise ValueError(f"Unknown candidate method: {method}")


def build_pattern_candidate_rows(
    categories: dict[str, list[str]],
    min_fragment_atoms: int = 3,
    min_prevalence: int = 2,
    max_breadth: int = 3,
    min_enrichment: float = 1.25,
) -> list[CandidateRecord]:
    methods = ("fixed_smarts", "murcko", "brics", "hybrid")
    method_category_counts: dict[str, dict[str, Counter[str]]] = {
        method: {category: Counter() for category in categories} for method in methods
    }

    for category, smiles_list in categories.items():
        for smiles in smiles_list:
            for method in methods:
                method_category_counts[method][category].update(extract_candidate_units(smiles, method, min_fragment_atoms))

    records: list[CandidateRecord] = []
    for method in methods:
        category_counts = method_category_counts[method]
        units = sorted({unit for counts in category_counts.values() for unit in counts})
        for unit in units:
            breadth = sum(1 for counts in category_counts.values() if counts.get(unit, 0) > 0)
            if breadth > max_breadth:
                continue
            for category in sorted(categories):
                count = category_counts[category].get(unit, 0)
                total = len(categories[category])
                if count < min_prevalence or total <= 0:
                    continue
                prevalence = count / total
                other_prevalence = max(
                    (category_counts[other].get(unit, 0) / len(categories[other]) for other in categories if other != category and len(categories[other]) > 0),
                    default=0.0,
                )
                enrichment = prevalence / max(other_prevalence, 1.0 / max(total, 1))
                if enrichment < min_enrichment:
                    continue
                records.append(CandidateRecord(method, category, unit, count, prevalence, total, breadth, enrichment))
    return sorted(records, key=lambda item: (item.method, item.category, -item.enrichment, -item.prevalence_count, item.unit))


def build_pattern_candidate_summary_rows(
    categories: dict[str, list[str]],
    records: list[CandidateRecord],
    min_fragment_atoms: int = 3,
    min_prevalence: int = 2,
    invalid_smiles_skipped: int = 0,
) -> list[dict[str, str]]:
    methods = ("fixed_smarts", "murcko", "brics", "hybrid")
    retained_by_method: dict[str, set[str]] = {
        method: {record.unit for record in records if record.method == method}
        for method in methods
    }
    occurrence_sets, breadth_by_method = _retained_unit_occurrence_sets(
        categories,
        retained_by_method,
        min_fragment_atoms,
        min_prevalence,
    )

    total_molecules = sum(len(values) for values in categories.values())
    rows: list[dict[str, str]] = []
    for method in methods:
        retained = retained_by_method[method]
        coverage_count = 0
        for smiles_list in categories.values():
            for smiles in smiles_list:
                if extract_candidate_units(smiles, method, min_fragment_atoms) & retained:
                    coverage_count += 1
        unit_count = len(retained)
        exclusive_count = sum(1 for unit in retained if breadth_by_method[method][unit] == 1)
        shared_count = sum(1 for unit in retained if breadth_by_method[method][unit] > 1)
        shared_three_plus_count = sum(1 for unit in retained if breadth_by_method[method][unit] >= 3)
        mean_category_breadth = (
            sum(breadth_by_method[method][unit] for unit in retained) / unit_count if unit_count else 0.0
        )
        rows.append(
            {
                "method": method,
                "unique_retained_units": str(unit_count),
                "molecule_coverage_count": str(coverage_count),
                "molecule_coverage_fraction": _fmt(coverage_count / total_molecules if total_molecules else 0.0),
                "exclusive_unit_fraction": _fmt(exclusive_count / unit_count if unit_count else 0.0),
                "shared_unit_fraction": _fmt(shared_count / unit_count if unit_count else 0.0),
                "shared_three_plus_unit_fraction": _fmt(shared_three_plus_count / unit_count if unit_count else 0.0),
                "mean_category_breadth": _fmt(mean_category_breadth),
                "mean_pairwise_jaccard": _fmt(_mean_pairwise_jaccard(list(occurrence_sets[method].values()))),
                "category_count": str(len(categories)),
                "molecule_count": str(total_molecules),
                "invalid_smiles_skipped": str(invalid_smiles_skipped),
                "claim_scope": "diversity_comparison_only_no_performance_claim",
            }
        )
    return rows


def _retained_unit_occurrence_sets(
    categories: dict[str, list[str]],
    retained_by_method: dict[str, set[str]],
    min_fragment_atoms: int,
    min_prevalence: int,
) -> tuple[dict[str, dict[str, set[str]]], dict[str, Counter[str]]]:
    occurrence_sets: dict[str, dict[str, set[str]]] = {
        method: {category: set() for category in categories}
        for method in retained_by_method
    }
    breadth_by_method: dict[str, Counter[str]] = {method: Counter() for method in retained_by_method}
    for method, retained_units in retained_by_method.items():
        if not retained_units:
            continue
        for category, smiles_list in categories.items():
            counts: Counter[str] = Counter()
            for smiles in smiles_list:
                counts.update(extract_candidate_units(smiles, method, min_fragment_atoms) & retained_units)
            occurrence_sets[method][category] = {
                unit for unit, count in counts.items() if count >= min_prevalence
            }
        for unit in retained_units:
            breadth_by_method[method][unit] = sum(
                1 for category_units in occurrence_sets[method].values() if unit in category_units
            )
    return occurrence_sets, breadth_by_method


def _mean_pairwise_jaccard(unit_sets: list[set[str]]) -> float:
    if len(unit_sets) < 2:
        return 0.0
    values: list[float] = []
    for left_index, left in enumerate(unit_sets):
        for right in unit_sets[left_index + 1:]:
            union = left | right
            values.append((len(left & right) / len(union)) if union else 0.0)
    return sum(values) / len(values) if values else 0.0


def write_pattern_candidates(
    output_csv: str | Path,
    summary_csv: str | Path | None = None,
    input_csv: str | Path | None = None,
    category_column: str = "category",
    smiles_column: str | None = None,
    use_final_rebuild_inputs: bool = False,
    limit_per_category: int | None = None,
    min_fragment_atoms: int = 3,
    min_prevalence: int = 2,
    max_breadth: int = 3,
    min_enrichment: float = 1.25,
) -> Path:
    if use_final_rebuild_inputs:
        categories, invalid_count = read_final_rebuild_inputs(limit_per_category)
    elif input_csv:
        categories, invalid_count = read_category_smiles_csv(input_csv, category_column, smiles_column)
    else:
        raise ValueError("Provide --input-csv or --use-final-rebuild-inputs.")
    if not categories:
        raise ValueError("No valid category SMILES were found in the selected input.")
    records = build_pattern_candidate_rows(categories, min_fragment_atoms, min_prevalence, max_breadth, min_enrichment)
    rows = [
        {
            "method": record.method,
            "category": record.category,
            "unit": record.unit,
            "prevalence_count": str(record.prevalence_count),
            "prevalence_fraction": _fmt(record.prevalence_fraction),
            "category_total": str(record.category_total),
            "breadth": str(record.breadth),
            "enrichment": _fmt(record.enrichment),
            "invalid_smiles_skipped": str(invalid_count),
        }
        for record in records
    ]
    fields = [
        "method",
        "category",
        "unit",
        "prevalence_count",
        "prevalence_fraction",
        "category_total",
        "breadth",
        "enrichment",
        "invalid_smiles_skipped",
    ]
    output_path = Path(output_csv)
    _write_csv(output_path, rows, fields)
    if summary_csv is not None:
        summary_rows = build_pattern_candidate_summary_rows(categories, records, min_fragment_atoms, min_prevalence, invalid_count)
        summary_fields = [
            "method",
            "unique_retained_units",
            "molecule_coverage_count",
            "molecule_coverage_fraction",
            "exclusive_unit_fraction",
            "shared_unit_fraction",
            "shared_three_plus_unit_fraction",
            "mean_category_breadth",
            "mean_pairwise_jaccard",
            "category_count",
            "molecule_count",
            "invalid_smiles_skipped",
            "claim_scope",
        ]
        _write_csv(Path(summary_csv), summary_rows, summary_fields)
    return output_path


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Deterministic public-model audit utilities.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    audit = subparsers.add_parser("audit", help="Write seven-probe public-model audit CSVs.")
    audit.add_argument("--panel", default=str(DEFAULT_PROBE_PANEL))
    audit.add_argument("--scores-out", default=str(DEFAULT_AUDIT_SCORES))
    audit.add_argument("--summary-out", default=str(DEFAULT_AUDIT_SUMMARY))

    overlap = subparsers.add_parser("pattern-overlap", help="Write public-model fixed-pattern overlap report.")
    overlap.add_argument("--csv-out", default=str(DEFAULT_PATTERN_OVERLAP))
    overlap.add_argument("--markdown-out", default=str(DEFAULT_PATTERN_AUDIT))

    candidates = subparsers.add_parser("pattern-candidates", help="Compare fixed SMARTS, Murcko, BRICS, and hybrid pattern candidates.")
    candidates.add_argument("--output-csv", required=True)
    candidates.add_argument("--summary-csv", help="Defaults to <output stem>_summary.csv.")
    candidates.add_argument("--input-csv")
    candidates.add_argument("--category-column", default="category")
    candidates.add_argument("--smiles-column")
    candidates.add_argument("--use-final-rebuild-inputs", action="store_true")
    candidates.add_argument("--limit-per-category", type=int)
    candidates.add_argument("--min-fragment-atoms", type=int, default=3)
    candidates.add_argument("--min-prevalence", type=int, default=2)
    candidates.add_argument("--max-breadth", type=int, default=3)
    candidates.add_argument("--min-enrichment", type=float, default=1.25)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    if args.command == "audit":
        write_default_probe_panel(args.panel)
        scores, summary = run_audit(args.panel, args.scores_out, args.summary_out)
        print(f"Wrote {scores}")
        print(f"Wrote {summary}")
        return
    if args.command == "pattern-overlap":
        csv_output, markdown_output = write_pattern_overlap(args.csv_out, args.markdown_out)
        print(f"Wrote {csv_output}")
        print(f"Wrote {markdown_output}")
        return
    if args.command == "pattern-candidates":
        summary_csv = args.summary_csv or str(Path(args.output_csv).with_name(f"{Path(args.output_csv).stem}_summary.csv"))
        output = write_pattern_candidates(
            output_csv=args.output_csv,
            summary_csv=summary_csv,
            input_csv=args.input_csv,
            category_column=args.category_column,
            smiles_column=args.smiles_column,
            use_final_rebuild_inputs=args.use_final_rebuild_inputs,
            limit_per_category=args.limit_per_category,
            min_fragment_atoms=args.min_fragment_atoms,
            min_prevalence=args.min_prevalence,
            max_breadth=args.max_breadth,
            min_enrichment=args.min_enrichment,
        )
        print(f"Wrote {output}")
        print(f"Wrote {summary_csv}")
        return
    raise AssertionError(args.command)


if __name__ == "__main__":
    main()
