"""Leakage-resistant loop for structural-pattern discovery and validation.

The loop is deliberately research-only.  It writes staged model candidates and
validation artifacts under an output directory, but never changes bundled models.
Candidate patterns are mined from training positives, selected against training
background categories, tuned on a validation split, and evaluated once with the
frozen threshold on globally scaffold-disjoint test molecules.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import BRICS, Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

try:
    from app.structural_pattern_library import CHOI_CANDIDATE_PATTERNS
except ModuleNotFoundError:  # Preserve direct execution from the repository.
    from structural_pattern_library import CHOI_CANDIDATE_PATTERNS


RDLogger.DisableLog("rdApp.*")

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
DEFAULT_INPUT_DIR = APP_DIR / "output" / "final_category_rebuild" / "inputs"
DEFAULT_OUTPUT_DIR = ROOT_DIR / "results" / "structural_pattern_validation"

PRODUCT_USE_CATEGORIES = (
    "animal_drugs",
    "cosmetics",
    "flavoring_agents",
    "food_additives",
    "food_contact_substances",
    "fragrances",
    "human_drugs",
    "pesticides",
    "solvents",
    "surfactants",
)
STRUCTURAL_METHODS = ("fixed_smarts", "murcko", "brics", "hybrid")
CANDIDATE_METHODS = ("murcko", "brics", "hybrid")
EVALUATION_REGIMES = ("all_other", "related_hard", "property_matched")
SCHEMA_VERSION = "structural-pattern-validation-v1"

PROPERTY_NAMES = (
    "MW",
    "logP",
    "HBD",
    "HBA",
    "TPSA",
    "RotBonds",
    "FCsp3",
    "AromaticRings",
)
MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@dataclass(frozen=True)
class ValidationConfig:
    input_dir: Path = DEFAULT_INPUT_DIR
    output_dir: Path = DEFAULT_OUTPUT_DIR
    seeds: tuple[int, ...] = (11, 23, 37)
    categories: tuple[str, ...] = PRODUCT_USE_CATEGORIES
    limit_per_category: int | None = None
    train_fraction: float = 0.60
    validation_fraction: float = 0.20
    min_split_size: int = 2
    min_fragment_atoms: int = 3
    max_fragment_atoms: int = 18
    candidate_pool_size: int = 96
    max_patterns: int = 24
    min_positive_count: int = 3
    min_positive_prevalence: float = 0.01
    min_enrichment: float = 1.25
    min_specificity_gap: float = 0.01
    min_cross_category_prevalence: float = 0.01
    max_category_breadth: int = 6
    property_ks_threshold: float = 0.15
    property_lower_quantile: float = 0.10
    property_upper_quantile: float = 0.90
    weight_grid: tuple[float, ...] = (0.0, 0.25, 0.50, 0.75, 1.0)
    max_all_other_negatives: int = 2400
    hard_negative_limit: int = 600
    similarity_reference_limit: int = 96
    bootstrap_replicates: int = 300
    minimum_promotion_seeds: int = 3
    promotion_min_worst_auc: float = 0.80
    promotion_min_worst_balanced_accuracy: float = 0.75
    promotion_min_delta_auc: float = 0.01
    promotion_min_delta_balanced_accuracy: float = 0.01
    promotion_min_specificity_gap: float = 0.10
    promotion_max_off_target_rate: float = 0.60
    promotion_min_method_consistency: float = 2.0 / 3.0
    promotion_min_pattern_stability: float = 0.20
    promotion_require_positive_ci: bool = True

    def validate(self) -> None:
        if not self.seeds:
            raise ValueError("At least one seed is required.")
        if not self.categories:
            raise ValueError("At least one target category is required.")
        unknown = sorted(set(self.categories) - set(PRODUCT_USE_CATEGORIES))
        if unknown:
            raise ValueError(f"Unknown product-use categories: {', '.join(unknown)}")
        if not 0.0 < self.train_fraction < 1.0:
            raise ValueError("train_fraction must be between 0 and 1.")
        if not 0.0 < self.validation_fraction < 1.0:
            raise ValueError("validation_fraction must be between 0 and 1.")
        if self.train_fraction + self.validation_fraction >= 1.0:
            raise ValueError("train_fraction + validation_fraction must be less than 1.")
        if self.min_fragment_atoms < 1 or self.max_fragment_atoms < self.min_fragment_atoms:
            raise ValueError("Invalid fragment atom limits.")
        if self.max_patterns < 1 or self.candidate_pool_size < self.max_patterns:
            raise ValueError("candidate_pool_size must be at least max_patterns.")
        if self.min_positive_count < 1 or self.min_split_size < 1:
            raise ValueError("Count thresholds must be positive.")
        if self.max_all_other_negatives < 2 or self.hard_negative_limit < 2:
            raise ValueError("Negative regime limits must be at least two.")
        if self.bootstrap_replicates < 0:
            raise ValueError("bootstrap_replicates cannot be negative.")
        if not 0.0 <= self.property_lower_quantile < self.property_upper_quantile <= 1.0:
            raise ValueError("Invalid property quantile interval.")
        if any(not 0.0 <= value <= 1.0 for value in self.weight_grid):
            raise ValueError("weight_grid values must be between 0 and 1.")
        validate_output_location(self.output_dir, self.input_dir)


@dataclass(frozen=True)
class MoleculeRecord:
    category: str
    smiles: str
    smiles_hash: str
    scaffold_key: str
    scaffold_hash: str
    split: str
    props: tuple[float, ...]
    mol: Chem.Mol = field(compare=False, repr=False)


@dataclass(frozen=True)
class CandidateSpec:
    pattern_id: str
    method: str
    origin: str
    source_unit: str
    smarts: str
    discovery_count: int


@dataclass(frozen=True)
class PatternAssessment:
    spec: CandidateSpec
    positive_count: int
    positive_rate: float
    negative_rate: float
    max_other_rate: float
    category_breadth: int
    enrichment: float
    specificity_gap: float
    quality: float
    selected: bool
    rejection_reason: str
    weight: float


@dataclass
class FrozenModel:
    category: str
    method: str
    selected_props: tuple[str, ...]
    ranges: dict[str, tuple[float, float]]
    patterns: dict[str, str]
    pattern_weights: dict[str, float]
    best_w: float
    threshold: float
    validation_selection_score: float
    validation_metrics: dict[str, dict[str, float]]
    compiled_patterns: dict[str, Chem.Mol] = field(default_factory=dict, repr=False)

    def compile(self) -> None:
        compiled: dict[str, Chem.Mol] = {}
        for pattern_id, smarts in self.patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                raise ValueError(f"Candidate SMARTS does not compile: {pattern_id}={smarts!r}")
            compiled[pattern_id] = pattern
        if set(compiled) != set(self.pattern_weights):
            raise ValueError("pattern_weights keys must exactly match selected_patterns keys.")
        self.compiled_patterns = compiled


def canonical_smiles(value: str) -> str:
    mol = Chem.MolFromSmiles(value.strip()) if value and value.strip() else None
    return Chem.MolToSmiles(mol, canonical=True) if mol is not None else ""


def stable_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def validate_output_location(output_dir: Path, input_dir: Path) -> None:
    output = output_dir.resolve()
    inputs = input_dir.resolve()
    repository = ROOT_DIR.resolve()
    allowed_repository_roots = (
        (ROOT_DIR / "results").resolve(),
        (APP_DIR / "output").resolve(),
    )
    if _is_within(output, inputs) or _is_within(inputs, output):
        raise ValueError(
            f"output_dir must not contain or be contained by the input directory: {output}"
        )
    if _is_within(output, repository) and not any(
        _is_within(output, allowed) for allowed in allowed_repository_roots
    ):
        raise ValueError(
            "Repository-internal output_dir must be under results/ or app/output/; "
            f"refusing protected source, deployed-model, or manuscript path: {output}"
        )


def deterministic_sample(values: Sequence[Any], limit: int | None, salt: str = "") -> list[Any]:
    items = list(values)
    if limit is None or len(items) <= limit:
        return sorted(items, key=lambda item: str(item))
    ranked = sorted((stable_hash(f"{salt}|{item}"), item) for item in items)
    return sorted((item for _digest, item in ranked[:limit]), key=lambda item: str(item))


def _resolve_smiles_column(fieldnames: Sequence[str]) -> str:
    normalized = {name.strip().lower(): name for name in fieldnames}
    for candidate in ("smiles", "canonical_smiles", "isomeric_smiles"):
        if candidate in normalized:
            return normalized[candidate]
    raise ValueError(f"No SMILES column found in: {', '.join(fieldnames)}")


def read_smiles_csv(path: Path, limit: int | None = None) -> tuple[list[str], int]:
    values: set[str] = set()
    invalid = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Input CSV has no header: {path}")
        smiles_column = _resolve_smiles_column(reader.fieldnames)
        for row in reader:
            value = canonical_smiles(str(row.get(smiles_column, "")))
            if value:
                values.add(value)
            else:
                invalid += 1
    return deterministic_sample(sorted(values), limit, path.name), invalid


def load_category_inputs(config: ValidationConfig) -> tuple[dict[str, list[str]], dict[str, Any]]:
    categories: dict[str, list[str]] = {}
    input_rows: list[dict[str, Any]] = []
    total_invalid = 0
    # Always load every product-use category: non-target categories define the
    # one-vs-constructed-background evaluation for each requested target.
    for category in PRODUCT_USE_CATEGORIES:
        path = config.input_dir / f"{category}__positive.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing final rebuild input: {path}")
        smiles, invalid = read_smiles_csv(path, config.limit_per_category)
        if not smiles:
            raise ValueError(f"No valid SMILES found in: {path}")
        categories[category] = smiles
        total_invalid += invalid
        input_rows.append(
            {
                "category": category,
                "path": str(path.resolve()),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "valid_unique_count": len(smiles),
                "invalid_count": invalid,
            }
        )
    return categories, {"files": input_rows, "invalid_count": total_invalid}


def molecule_properties(mol: Chem.Mol) -> tuple[float, ...]:
    return (
        float(Descriptors.MolWt(mol)),
        float(Crippen.MolLogP(mol)),
        float(rdMolDescriptors.CalcNumHBD(mol)),
        float(rdMolDescriptors.CalcNumHBA(mol)),
        float(rdMolDescriptors.CalcTPSA(mol)),
        float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        float(rdMolDescriptors.CalcFractionCSP3(mol)),
        float(rdMolDescriptors.CalcNumAromaticRings(mol)),
    )


def _clean_brics_fragment(
    fragment: str,
    parent_smiles: str,
    min_fragment_atoms: int,
    max_fragment_atoms: int,
) -> str:
    mol = Chem.MolFromSmiles(fragment)
    if mol is None:
        return ""
    editable = Chem.RWMol(mol)
    dummy_indices = sorted(
        (atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0),
        reverse=True,
    )
    for atom_index in dummy_indices:
        editable.RemoveAtom(atom_index)
    clean = editable.GetMol()
    try:
        Chem.SanitizeMol(clean)
    except Exception:
        return ""
    heavy_atoms = clean.GetNumHeavyAtoms()
    if heavy_atoms < min_fragment_atoms or heavy_atoms > max_fragment_atoms:
        return ""
    value = Chem.MolToSmiles(clean, canonical=True)
    if not value or value == parent_smiles:
        return ""
    return value


def brics_fragment_units(
    mol: Chem.Mol,
    min_fragment_atoms: int = 3,
    max_fragment_atoms: int = 18,
) -> set[str]:
    parent = Chem.MolToSmiles(mol, canonical=True)
    return {
        cleaned
        for fragment in BRICS.BRICSDecompose(mol)
        if (
            cleaned := _clean_brics_fragment(
                fragment,
                parent,
                min_fragment_atoms,
                max_fragment_atoms,
            )
        )
    }


def murcko_fragment_units(
    mol: Chem.Mol,
    min_fragment_atoms: int = 3,
    max_fragment_atoms: int = 18,
) -> set[str]:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return set()
    heavy_atoms = scaffold.GetNumHeavyAtoms()
    if heavy_atoms < min_fragment_atoms or heavy_atoms > max_fragment_atoms:
        return set()
    value = Chem.MolToSmiles(scaffold, canonical=True)
    return {value} if value else set()


def scaffold_group_key(mol: Chem.Mol) -> str:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is not None and scaffold.GetNumHeavyAtoms() > 0:
        return f"murcko:{Chem.MolToSmiles(scaffold, canonical=True)}"

    # Bemis-Murcko is empty for acyclic molecules.  Genericize the full
    # molecular topology so atom/bond substitutions with the same topology stay
    # grouped without running fragment discovery before the split is frozen.
    try:
        generic = MurckoScaffold.MakeScaffoldGeneric(mol)
        generic_value = Chem.MolToSmiles(generic, canonical=True)
    except Exception:
        generic_value = ""
    if generic_value:
        return f"acyclic_generic:{generic_value}"
    return f"acyclic_molecule:{Chem.MolToSmiles(mol, canonical=True)}"


def partition_for_scaffold(
    scaffold_key: str,
    seed: int,
    train_fraction: float,
    validation_fraction: float,
) -> str:
    digest = stable_hash(f"{seed}|{scaffold_key}")
    fraction = int(digest[:16], 16) / float(16**16)
    if fraction < train_fraction:
        return "train"
    if fraction < train_fraction + validation_fraction:
        return "validation"
    return "test"


def build_split_records(
    categories: dict[str, Sequence[str]],
    seed: int,
    train_fraction: float = 0.60,
    validation_fraction: float = 0.20,
    molecule_cache: dict[str, tuple[Chem.Mol, str, tuple[float, ...]]] | None = None,
) -> dict[str, dict[str, list[MoleculeRecord]]]:
    records: dict[str, dict[str, list[MoleculeRecord]]] = {
        category: {"train": [], "validation": [], "test": []}
        for category in sorted(categories)
    }
    cache = molecule_cache if molecule_cache is not None else {}
    for category in sorted(categories):
        for smiles in sorted(set(categories[category])):
            cached = cache.get(smiles)
            if cached is None:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                scaffold_key = scaffold_group_key(mol)
                cached = (mol, scaffold_key, molecule_properties(mol))
                cache[smiles] = cached
            mol, scaffold_key, props = cached
            split = partition_for_scaffold(scaffold_key, seed, train_fraction, validation_fraction)
            record = MoleculeRecord(
                category=category,
                smiles=smiles,
                smiles_hash=stable_hash(smiles),
                scaffold_key=scaffold_key,
                scaffold_hash=stable_hash(scaffold_key),
                split=split,
                props=props,
                mol=mol,
            )
            records[category][split].append(record)

    for split_map in records.values():
        for split in split_map:
            split_map[split].sort(key=lambda row: row.smiles)
    assert_global_scaffold_disjoint(records)
    return records


def prepare_molecule_cache(
    categories: dict[str, Sequence[str]],
) -> dict[str, tuple[Chem.Mol, str, tuple[float, ...]]]:
    cache: dict[str, tuple[Chem.Mol, str, tuple[float, ...]]] = {}
    for smiles in sorted({value for values in categories.values() for value in values}):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        cache[smiles] = (mol, scaffold_group_key(mol), molecule_properties(mol))
    return cache


def assert_global_scaffold_disjoint(records: dict[str, dict[str, list[MoleculeRecord]]]) -> None:
    scaffold_splits: dict[str, set[str]] = defaultdict(set)
    molecule_splits: dict[str, set[str]] = defaultdict(set)
    for split_map in records.values():
        for split, rows in split_map.items():
            for row in rows:
                scaffold_splits[row.scaffold_key].add(split)
                molecule_splits[row.smiles].add(split)
    leaking_scaffolds = sorted(key for key, splits in scaffold_splits.items() if len(splits) > 1)
    leaking_molecules = sorted(key for key, splits in molecule_splits.items() if len(splits) > 1)
    if leaking_scaffolds or leaking_molecules:
        raise AssertionError(
            f"Global split leakage: {len(leaking_scaffolds)} scaffolds, "
            f"{len(leaking_molecules)} canonical molecules"
        )


def validate_split_sizes(
    records: dict[str, dict[str, list[MoleculeRecord]]],
    target_categories: Iterable[str],
    minimum: int,
) -> None:
    failures: list[str] = []
    for category in target_categories:
        for split in ("train", "validation", "test"):
            count = len(records[category][split])
            if count < minimum:
                failures.append(f"{category}/{split}={count}")
    if failures:
        raise ValueError(
            "Scaffold split is too small for validation (change seed or increase data): "
            + ", ".join(failures)
        )


def unique_records(rows: Iterable[MoleculeRecord], excluded_smiles: set[str] | None = None) -> list[MoleculeRecord]:
    excluded = excluded_smiles or set()
    by_smiles: dict[str, MoleculeRecord] = {}
    for row in rows:
        if row.smiles not in excluded:
            by_smiles.setdefault(row.smiles, row)
    return [by_smiles[key] for key in sorted(by_smiles)]


def background_by_category(
    records: dict[str, dict[str, list[MoleculeRecord]]],
    target_category: str,
    split: str,
) -> dict[str, list[MoleculeRecord]]:
    target_members = {
        row.smiles
        for target_split in records[target_category].values()
        for row in target_split
    }
    return {
        category: unique_records(records[category][split], target_members)
        for category in sorted(records)
        if category != target_category
    }


def flatten_background(by_category: dict[str, list[MoleculeRecord]]) -> list[MoleculeRecord]:
    return unique_records(row for rows in by_category.values() for row in rows)


def sample_background_by_category(
    background: dict[str, list[MoleculeRecord]],
    total_limit: int,
    salt: str,
) -> dict[str, list[MoleculeRecord]]:
    if not background:
        return {}
    per_category = max(1, math.ceil(total_limit / len(background)))
    return {
        category: _sample_records(rows, per_category, f"{salt}|{category}")
        for category, rows in sorted(background.items())
    }


def _sample_records(rows: Sequence[MoleculeRecord], limit: int | None, salt: str) -> list[MoleculeRecord]:
    if limit is None or len(rows) <= limit:
        return sorted(rows, key=lambda row: row.smiles)
    ranked = sorted((stable_hash(f"{salt}|{row.smiles}"), row.smiles, row) for row in rows)
    return sorted((row for _digest, _smiles, row in ranked[:limit]), key=lambda row: row.smiles)


def _unit_to_smarts(unit: str) -> str:
    mol = Chem.MolFromSmiles(unit)
    if mol is None or any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        return ""
    smarts = Chem.MolToSmarts(mol)
    return smarts if smarts and Chem.MolFromSmarts(smarts) is not None else ""


def _canonical_query_smarts(smarts: str) -> str:
    query = Chem.MolFromSmarts(smarts)
    return Chem.MolToSmarts(query) if query is not None else ""


def discover_candidate_specs(
    train_positive: Sequence[MoleculeRecord],
    method: str,
    config: ValidationConfig,
) -> list[CandidateSpec]:
    if method not in STRUCTURAL_METHODS:
        raise ValueError(f"Unknown structural method: {method}")

    raw: list[tuple[str, str, str, int]] = []
    if method in ("fixed_smarts", "hybrid"):
        for name, smarts in sorted(CHOI_CANDIDATE_PATTERNS.items()):
            query = Chem.MolFromSmarts(smarts)
            if query is None:
                continue
            count = sum(row.mol.HasSubstructMatch(query) for row in train_positive)
            raw.append(("fixed_smarts", name, _canonical_query_smarts(smarts), int(count)))

    generated_methods = (
        ("murcko", murcko_fragment_units),
        ("brics", brics_fragment_units),
    )
    for origin, extractor in generated_methods:
        if method not in (origin, "hybrid"):
            continue
        counts: Counter[str] = Counter()
        for row in train_positive:
            counts.update(
                extractor(
                    row.mol,
                    min_fragment_atoms=config.min_fragment_atoms,
                    max_fragment_atoms=config.max_fragment_atoms,
                )
            )
        eligible = [
            (unit, count)
            for unit, count in counts.items()
            if count >= config.min_positive_count
            and count / max(len(train_positive), 1) >= config.min_positive_prevalence
        ]
        eligible.sort(key=lambda item: (-item[1], item[0]))
        for unit, count in eligible[: config.candidate_pool_size]:
            smarts = _unit_to_smarts(unit)
            if smarts:
                raw.append((origin, unit, smarts, int(count)))

    # Collapse equivalent queries before assigning stable IDs.  In hybrid mode,
    # fixed SMARTS win ties so familiar motifs keep their human-readable source.
    origin_order = {"fixed_smarts": 0, "brics": 1, "murcko": 2}
    by_smarts: dict[str, tuple[str, str, int]] = {}
    for origin, source_unit, smarts, count in sorted(
        raw,
        key=lambda item: (origin_order[item[0]], -item[3], item[1]),
    ):
        previous = by_smarts.get(smarts)
        if previous is None or count > previous[2]:
            by_smarts[smarts] = (origin, source_unit, count)

    specs: list[CandidateSpec] = []
    for smarts, (origin, source_unit, count) in sorted(by_smarts.items()):
        digest = stable_hash(smarts)[:12]
        readable = source_unit if origin == "fixed_smarts" else digest
        pattern_id = f"{origin}__{readable}"
        specs.append(
            CandidateSpec(
                pattern_id=pattern_id,
                method=method,
                origin=origin,
                source_unit=source_unit,
                smarts=smarts,
                discovery_count=count,
            )
        )
    return sorted(specs, key=lambda item: item.pattern_id)


def assess_and_select_patterns(
    specs: Sequence[CandidateSpec],
    train_positive: Sequence[MoleculeRecord],
    train_background: dict[str, list[MoleculeRecord]],
    config: ValidationConfig,
) -> tuple[list[PatternAssessment], dict[str, str], dict[str, float]]:
    negative = flatten_background(train_background)
    assessments: list[PatternAssessment] = []
    for spec in specs:
        query = Chem.MolFromSmarts(spec.smarts)
        if query is None:
            continue
        positive_count = int(sum(row.mol.HasSubstructMatch(query) for row in train_positive))
        positive_rate = positive_count / max(len(train_positive), 1)
        negative_count = int(sum(row.mol.HasSubstructMatch(query) for row in negative))
        negative_rate = negative_count / max(len(negative), 1)
        other_rates = {
            category: sum(row.mol.HasSubstructMatch(query) for row in rows) / max(len(rows), 1)
            for category, rows in train_background.items()
        }
        max_other_rate = max(other_rates.values(), default=0.0)
        breadth = sum(
            rate >= config.min_cross_category_prevalence
            for rate in other_rates.values()
        )
        smoothed_positive = (positive_count + 0.5) / (len(train_positive) + 1.0)
        smoothed_negative = (negative_count + 0.5) / (len(negative) + 1.0)
        enrichment = smoothed_positive / max(smoothed_negative, 1e-12)
        specificity_gap = positive_rate - max_other_rate
        quality = max(specificity_gap, 0.0) * math.log2(1.0 + enrichment) * math.sqrt(max(positive_rate, 0.0))

        reasons: list[str] = []
        if positive_count < config.min_positive_count:
            reasons.append("positive_count")
        if positive_rate < config.min_positive_prevalence:
            reasons.append("positive_prevalence")
        if enrichment < config.min_enrichment:
            reasons.append("enrichment")
        if specificity_gap < config.min_specificity_gap:
            reasons.append("specificity_gap")
        if breadth > config.max_category_breadth:
            reasons.append("category_breadth")
        if quality <= 0.0:
            reasons.append("nonpositive_quality")

        assessments.append(
            PatternAssessment(
                spec=spec,
                positive_count=positive_count,
                positive_rate=float(positive_rate),
                negative_rate=float(negative_rate),
                max_other_rate=float(max_other_rate),
                category_breadth=int(breadth),
                enrichment=float(enrichment),
                specificity_gap=float(specificity_gap),
                quality=float(quality),
                selected=not reasons,
                rejection_reason=";".join(sorted(set(reasons))),
                weight=float(quality),
            )
        )

    eligible = sorted(
        (item for item in assessments if item.selected),
        key=lambda item: (-item.quality, -item.positive_rate, item.spec.pattern_id),
    )
    selected_ids = {item.spec.pattern_id for item in eligible[: config.max_patterns]}
    final_assessments: list[PatternAssessment] = []
    for item in assessments:
        if item.selected and item.spec.pattern_id not in selected_ids:
            item = PatternAssessment(
                spec=item.spec,
                positive_count=item.positive_count,
                positive_rate=item.positive_rate,
                negative_rate=item.negative_rate,
                max_other_rate=item.max_other_rate,
                category_breadth=item.category_breadth,
                enrichment=item.enrichment,
                specificity_gap=item.specificity_gap,
                quality=item.quality,
                selected=False,
                rejection_reason="max_patterns",
                weight=0.0,
            )
        final_assessments.append(item)

    selected = [item for item in final_assessments if item.selected]
    patterns = {item.spec.pattern_id: item.spec.smarts for item in selected}
    weights = {item.spec.pattern_id: item.weight for item in selected}
    return sorted(final_assessments, key=lambda item: item.spec.pattern_id), patterns, weights


def ks_statistic(left: Sequence[float], right: Sequence[float]) -> float:
    a = np.sort(np.asarray(left, dtype=float))
    b = np.sort(np.asarray(right, dtype=float))
    if len(a) == 0 or len(b) == 0:
        return 0.0
    values = np.sort(np.unique(np.r_[a, b]))
    cdf_a = np.searchsorted(a, values, side="right") / len(a)
    cdf_b = np.searchsorted(b, values, side="right") / len(b)
    return float(np.max(np.abs(cdf_a - cdf_b)))


def fit_property_component(
    train_positive: Sequence[MoleculeRecord],
    train_negative: Sequence[MoleculeRecord],
    config: ValidationConfig,
) -> tuple[tuple[str, ...], dict[str, tuple[float, float]], dict[str, float]]:
    pos_values = np.asarray([row.props for row in train_positive], dtype=float)
    neg_values = np.asarray([row.props for row in train_negative], dtype=float)
    if pos_values.size == 0 or neg_values.size == 0:
        return (), {}, {}
    statistics = {
        name: ks_statistic(pos_values[:, index], neg_values[:, index])
        for index, name in enumerate(PROPERTY_NAMES)
    }
    selected = tuple(
        name for name in PROPERTY_NAMES if statistics[name] >= config.property_ks_threshold
    )
    if not selected and statistics:
        selected = (max(PROPERTY_NAMES, key=lambda name: (statistics[name], name)),)
    ranges = {
        name: (
            float(np.quantile(pos_values[:, PROPERTY_NAMES.index(name)], config.property_lower_quantile)),
            float(np.quantile(pos_values[:, PROPERTY_NAMES.index(name)], config.property_upper_quantile)),
        )
        for name in selected
    }
    return selected, ranges, statistics


def score_components(
    rows: Sequence[MoleculeRecord],
    model: FrozenModel,
) -> tuple[np.ndarray, np.ndarray]:
    if model.patterns and not model.compiled_patterns:
        model.compile()
    prop_indices = {name: PROPERTY_NAMES.index(name) for name in model.selected_props}
    total_weight = sum(model.pattern_weights.values()) or 1.0
    property_scores: list[float] = []
    structure_scores: list[float] = []
    for row in rows:
        if model.selected_props:
            property_score = sum(
                model.ranges[name][0] <= row.props[prop_indices[name]] <= model.ranges[name][1]
                for name in model.selected_props
            ) / len(model.selected_props)
        else:
            property_score = 0.0
        structure_score = sum(
            model.pattern_weights[pattern_id]
            for pattern_id, query in model.compiled_patterns.items()
            if row.mol.HasSubstructMatch(query)
        ) / total_weight
        property_scores.append(float(property_score))
        structure_scores.append(float(structure_score))
    return np.asarray(property_scores), np.asarray(structure_scores)


def score_records(rows: Sequence[MoleculeRecord], model: FrozenModel) -> np.ndarray:
    property_scores, structure_scores = score_components(rows, model)
    return model.best_w * property_scores + (1.0 - model.best_w) * structure_scores


def related_hard_negatives(
    candidates: Sequence[MoleculeRecord],
    train_positive: Sequence[MoleculeRecord],
    limit: int,
    reference_limit: int,
    salt: str,
) -> list[MoleculeRecord]:
    if not candidates or not train_positive:
        return []
    references = _sample_records(train_positive, reference_limit, f"{salt}|references")
    reference_fps = [MORGAN_GENERATOR.GetFingerprint(row.mol) for row in references]
    ranked: list[tuple[float, str, MoleculeRecord]] = []
    for row in candidates:
        fingerprint = MORGAN_GENERATOR.GetFingerprint(row.mol)
        maximum = max(DataStructs.BulkTanimotoSimilarity(fingerprint, reference_fps), default=0.0)
        ranked.append((-float(maximum), row.smiles, row))
    ranked.sort(key=lambda item: (item[0], item[1]))
    return [row for _negative_similarity, _smiles, row in ranked[:limit]]


def property_matched_negatives(
    candidates: Sequence[MoleculeRecord],
    train_positive: Sequence[MoleculeRecord],
    limit: int,
) -> list[MoleculeRecord]:
    if not candidates or not train_positive:
        return []
    positive = np.asarray([row.props for row in train_positive], dtype=float)
    center = np.median(positive, axis=0)
    scale = np.quantile(positive, 0.75, axis=0) - np.quantile(positive, 0.25, axis=0)
    scale = np.where(scale < 1e-8, np.std(positive, axis=0), scale)
    scale = np.where(scale < 1e-8, 1.0, scale)
    ranked = sorted(
        (
            float(np.linalg.norm((np.asarray(row.props) - center) / scale)),
            row.smiles,
            row,
        )
        for row in candidates
    )
    return [row for _distance, _smiles, row in ranked[:limit]]


def build_negative_regimes(
    background: dict[str, list[MoleculeRecord]],
    train_positive: Sequence[MoleculeRecord],
    config: ValidationConfig,
    salt: str,
) -> dict[str, list[MoleculeRecord]]:
    all_candidates = flatten_background(background)
    all_other = _sample_records(
        all_candidates,
        config.max_all_other_negatives,
        f"{salt}|all_other",
    )
    related_hard = related_hard_negatives(
        all_candidates,
        train_positive,
        min(config.hard_negative_limit, len(all_candidates)),
        config.similarity_reference_limit,
        f"{salt}|related_hard",
    )
    property_matched = property_matched_negatives(
        all_candidates,
        train_positive,
        min(config.hard_negative_limit, len(all_candidates)),
    )
    return {
        "all_other": all_other,
        "related_hard": related_hard,
        "property_matched": property_matched,
    }


def roc_auc(pos_scores: Sequence[float], neg_scores: Sequence[float]) -> float:
    positive = np.asarray(pos_scores, dtype=float)
    negative = np.asarray(neg_scores, dtype=float)
    if len(positive) == 0 or len(negative) == 0:
        return 0.0
    scores = np.r_[positive, negative]
    labels = np.r_[np.ones(len(positive), dtype=int), np.zeros(len(negative), dtype=int)]
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    index = 0
    while index < len(order):
        end = index + 1
        while end < len(order) and scores[order[end]] == scores[order[index]]:
            end += 1
        ranks[order[index:end]] = 0.5 * (index + 1 + end)
        index = end
    rank_sum = float(np.sum(ranks[labels == 1]))
    return float((rank_sum - len(positive) * (len(positive) + 1) / 2.0) / (len(positive) * len(negative)))


def average_precision(pos_scores: Sequence[float], neg_scores: Sequence[float]) -> float:
    positive = np.asarray(pos_scores, dtype=float)
    negative = np.asarray(neg_scores, dtype=float)
    if len(positive) == 0:
        return 0.0
    grouped: dict[float, list[int]] = defaultdict(lambda: [0, 0])
    for score in positive:
        grouped[float(score)][0] += 1
    for score in negative:
        grouped[float(score)][1] += 1
    true_positive = 0
    false_positive = 0
    previous_recall = 0.0
    area = 0.0
    for score in sorted(grouped, reverse=True):
        positives_at_score, negatives_at_score = grouped[score]
        true_positive += positives_at_score
        false_positive += negatives_at_score
        recall = true_positive / len(positive)
        precision = true_positive / max(true_positive + false_positive, 1)
        area += (recall - previous_recall) * precision
        previous_recall = recall
    return float(area)


def binary_metrics(
    pos_scores: Sequence[float],
    neg_scores: Sequence[float],
    threshold: float,
) -> dict[str, float]:
    positive = np.asarray(pos_scores, dtype=float)
    negative = np.asarray(neg_scores, dtype=float)
    if len(positive) == 0 or len(negative) == 0:
        return {
            "auc": 0.0,
            "average_precision": 0.0,
            "balanced_accuracy": 0.0,
            "sensitivity": 0.0,
            "specificity": 0.0,
            "mcc": 0.0,
            "ks": 0.0,
            "threshold": float(threshold),
            "positive_count": float(len(positive)),
            "negative_count": float(len(negative)),
        }
    true_positive = int(np.sum(positive >= threshold))
    false_negative = len(positive) - true_positive
    false_positive = int(np.sum(negative >= threshold))
    true_negative = len(negative) - false_positive
    sensitivity = true_positive / len(positive)
    specificity = true_negative / len(negative)
    denominator = math.sqrt(
        max(
            (true_positive + false_positive)
            * (true_positive + false_negative)
            * (true_negative + false_positive)
            * (true_negative + false_negative),
            0,
        )
    )
    mcc = (
        (true_positive * true_negative - false_positive * false_negative) / denominator
        if denominator > 0
        else 0.0
    )
    return {
        "auc": roc_auc(positive, negative),
        "average_precision": average_precision(positive, negative),
        "balanced_accuracy": float(0.5 * (sensitivity + specificity)),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "mcc": float(mcc),
        "ks": ks_statistic(positive, negative),
        "threshold": float(threshold),
        "positive_count": float(len(positive)),
        "negative_count": float(len(negative)),
    }


def choose_threshold(pos_scores: Sequence[float], neg_scores: Sequence[float]) -> tuple[float, dict[str, float]]:
    scores = np.r_[np.asarray(pos_scores, dtype=float), np.asarray(neg_scores, dtype=float)]
    if len(scores) == 0:
        return 0.5, binary_metrics([], [], 0.5)
    candidates = sorted(set(float(value) for value in scores))
    candidates.append(float(np.nextafter(max(candidates), math.inf)))
    best_threshold = candidates[0]
    best_metrics = binary_metrics(pos_scores, neg_scores, best_threshold)
    for threshold in candidates[1:]:
        metrics = binary_metrics(pos_scores, neg_scores, threshold)
        candidate_key = (
            metrics["balanced_accuracy"],
            metrics["mcc"],
            -abs(threshold - 0.5),
            -threshold,
        )
        best_key = (
            best_metrics["balanced_accuracy"],
            best_metrics["mcc"],
            -abs(best_threshold - 0.5),
            -best_threshold,
        )
        if candidate_key > best_key:
            best_threshold = threshold
            best_metrics = metrics
    return float(best_threshold), best_metrics


def cross_category_rates(
    model: FrozenModel,
    category_rows: dict[str, Sequence[MoleculeRecord]],
) -> dict[str, dict[str, float]]:
    rates: dict[str, dict[str, float]] = {}
    for category in sorted(category_rows):
        scores = score_records(category_rows[category], model)
        rates[category] = {
            "positive_rate": float(np.mean(scores >= model.threshold)) if len(scores) else 0.0,
            "mean_score": float(np.mean(scores)) if len(scores) else 0.0,
            "count": float(len(scores)),
        }
    return rates


def tune_frozen_model(
    category: str,
    method: str,
    selected_props: tuple[str, ...],
    ranges: dict[str, tuple[float, float]],
    patterns: dict[str, str],
    pattern_weights: dict[str, float],
    validation_positive: Sequence[MoleculeRecord],
    validation_negative_regimes: dict[str, list[MoleculeRecord]],
    validation_categories: dict[str, Sequence[MoleculeRecord]],
    weight_grid: Sequence[float],
) -> FrozenModel:
    allowed_weights = tuple(weight_grid)
    if method == "property_only":
        allowed_weights = (1.0,)
    elif method == "hybrid_structure_only":
        allowed_weights = (0.0,)
    elif not patterns:
        allowed_weights = (1.0,)

    best_model: FrozenModel | None = None
    for weight in allowed_weights:
        candidate = FrozenModel(
            category=category,
            method=method,
            selected_props=selected_props,
            ranges=ranges,
            patterns=dict(patterns),
            pattern_weights=dict(pattern_weights),
            best_w=float(weight),
            threshold=0.5,
            validation_selection_score=-math.inf,
            validation_metrics={},
        )
        candidate.compile()
        positive_scores = score_records(validation_positive, candidate)
        all_other_scores = score_records(validation_negative_regimes["all_other"], candidate)
        threshold, _metrics = choose_threshold(positive_scores, all_other_scores)
        candidate.threshold = threshold

        regime_metrics: dict[str, dict[str, float]] = {}
        for regime in EVALUATION_REGIMES:
            negative_scores = score_records(validation_negative_regimes[regime], candidate)
            regime_metrics[regime] = binary_metrics(positive_scores, negative_scores, threshold)

        rates = cross_category_rates(candidate, validation_categories)
        target_rate = rates.get(category, {}).get("positive_rate", 0.0)
        off_target_rates = [
            row["positive_rate"] for other, row in rates.items() if other != category
        ]
        max_off_target = max(off_target_rates, default=0.0)
        specificity_gap = target_rate - max_off_target
        worst_auc = min(row["auc"] for row in regime_metrics.values())
        worst_balanced_accuracy = min(row["balanced_accuracy"] for row in regime_metrics.values())
        selection_score = (
            0.40 * worst_balanced_accuracy
            + 0.35 * worst_auc
            + 0.25 * ((specificity_gap + 1.0) / 2.0)
        )
        candidate.validation_selection_score = float(selection_score)
        candidate.validation_metrics = regime_metrics
        candidate.validation_metrics["cross_category"] = {
            "target_positive_rate": float(target_rate),
            "max_off_target_rate": float(max_off_target),
            "specificity_gap": float(specificity_gap),
        }

        candidate_key = (
            candidate.validation_selection_score,
            worst_balanced_accuracy,
            worst_auc,
            -abs(weight - 0.5),
        )
        if best_model is None:
            best_model = candidate
        else:
            best_worst_auc = min(
                row["auc"]
                for name, row in best_model.validation_metrics.items()
                if name in EVALUATION_REGIMES
            )
            best_worst_ba = min(
                row["balanced_accuracy"]
                for name, row in best_model.validation_metrics.items()
                if name in EVALUATION_REGIMES
            )
            best_key = (
                best_model.validation_selection_score,
                best_worst_ba,
                best_worst_auc,
                -abs(best_model.best_w - 0.5),
            )
            if candidate_key > best_key:
                best_model = candidate
    if best_model is None:
        raise AssertionError("No validation model was evaluated.")
    return best_model


def evaluate_frozen_model(
    model: FrozenModel,
    positive: Sequence[MoleculeRecord],
    negative_regimes: dict[str, list[MoleculeRecord]],
) -> tuple[dict[str, dict[str, float]], dict[str, np.ndarray]]:
    positive_scores = score_records(positive, model)
    metrics: dict[str, dict[str, float]] = {}
    scores: dict[str, np.ndarray] = {"positive": positive_scores}
    for regime in EVALUATION_REGIMES:
        negative_scores = score_records(negative_regimes[regime], model)
        metrics[regime] = binary_metrics(positive_scores, negative_scores, model.threshold)
        metrics[regime]["threshold_source_validation"] = 1.0
        scores[regime] = negative_scores
    return metrics, scores


def paired_bootstrap_deltas(
    candidate_scores: dict[str, np.ndarray],
    candidate_threshold: float,
    baseline_scores: dict[str, np.ndarray],
    baseline_threshold: float,
    regime: str,
    replicates: int,
    seed: int,
) -> dict[str, float]:
    if replicates <= 0:
        return {
            "delta_auc_lower": float("nan"),
            "delta_auc_median": float("nan"),
            "delta_auc_upper": float("nan"),
            "delta_balanced_accuracy_lower": float("nan"),
            "delta_balanced_accuracy_median": float("nan"),
            "delta_balanced_accuracy_upper": float("nan"),
        }
    candidate_positive = candidate_scores["positive"]
    baseline_positive = baseline_scores["positive"]
    candidate_negative = candidate_scores[regime]
    baseline_negative = baseline_scores[regime]
    if len(candidate_positive) != len(baseline_positive) or len(candidate_negative) != len(baseline_negative):
        raise ValueError("Paired bootstrap requires aligned candidate and baseline scores.")

    generator = np.random.default_rng(seed)
    auc_deltas: list[float] = []
    balanced_accuracy_deltas: list[float] = []
    for _ in range(replicates):
        positive_indices = generator.integers(0, len(candidate_positive), len(candidate_positive))
        negative_indices = generator.integers(0, len(candidate_negative), len(candidate_negative))
        candidate_metrics = binary_metrics(
            candidate_positive[positive_indices],
            candidate_negative[negative_indices],
            candidate_threshold,
        )
        baseline_metrics = binary_metrics(
            baseline_positive[positive_indices],
            baseline_negative[negative_indices],
            baseline_threshold,
        )
        auc_deltas.append(candidate_metrics["auc"] - baseline_metrics["auc"])
        balanced_accuracy_deltas.append(
            candidate_metrics["balanced_accuracy"] - baseline_metrics["balanced_accuracy"]
        )

    def intervals(values: Sequence[float]) -> tuple[float, float, float]:
        return tuple(float(value) for value in np.quantile(values, (0.025, 0.50, 0.975)))  # type: ignore[return-value]

    auc_lower, auc_median, auc_upper = intervals(auc_deltas)
    ba_lower, ba_median, ba_upper = intervals(balanced_accuracy_deltas)
    return {
        "delta_auc_lower": auc_lower,
        "delta_auc_median": auc_median,
        "delta_auc_upper": auc_upper,
        "delta_balanced_accuracy_lower": ba_lower,
        "delta_balanced_accuracy_median": ba_median,
        "delta_balanced_accuracy_upper": ba_upper,
    }


def _model_to_json(
    model: FrozenModel,
    seed: int,
    test_metrics: dict[str, dict[str, float]],
    split_digest: str,
) -> dict[str, Any]:
    return {
        "model_id": f"research_{model.category}_{model.method}_seed_{seed}",
        "label": f"Research candidate: {model.category} ({model.method}, seed {seed})",
        "category": model.category,
        "model_type": "choi_auto",
        "threshold": model.threshold,
        "selected_props": list(model.selected_props),
        "ranges": {name: list(bounds) for name, bounds in sorted(model.ranges.items())},
        "selected_patterns": dict(sorted(model.patterns.items())),
        "pattern_weights": dict(sorted(model.pattern_weights.items())),
        "best_w": model.best_w,
        "metrics": test_metrics,
        "validation_metrics": model.validation_metrics,
        "optimization_method": "global_scaffold_train_validation_test",
        "threshold_source": "validation_all_other",
        "candidate_generation_source": "training_positives_only",
        "split_manifest_sha256": split_digest,
        "research_only": True,
        "automatic_release_allowed": False,
        "schema_version": SCHEMA_VERSION,
    }


def _split_rows(
    records: dict[str, dict[str, list[MoleculeRecord]]],
    seed: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for category in sorted(records):
        for split in ("train", "validation", "test"):
            for row in records[category][split]:
                rows.append(
                    {
                        "seed": seed,
                        "category": category,
                        "split": split,
                        "smiles_sha256": row.smiles_hash,
                        "scaffold_sha256": row.scaffold_hash,
                        "scaffold_type": row.scaffold_key.split(":", 1)[0],
                    }
                )
    return rows


def _split_digest(rows: Sequence[dict[str, Any]]) -> str:
    payload = json.dumps(list(rows), sort_keys=True, separators=(",", ":"))
    return stable_hash(payload)


def run_seed(
    categories: dict[str, list[str]],
    seed: int,
    config: ValidationConfig,
    molecule_cache: dict[str, tuple[Chem.Mol, str, tuple[float, ...]]] | None = None,
) -> dict[str, Any]:
    records = build_split_records(
        categories,
        seed,
        config.train_fraction,
        config.validation_fraction,
        molecule_cache,
    )
    validate_split_sizes(records, config.categories, config.min_split_size)
    split_rows = _split_rows(records, seed)
    split_digest = _split_digest(split_rows)

    metric_rows: list[dict[str, Any]] = []
    pattern_rows: list[dict[str, Any]] = []
    cross_rows: list[dict[str, Any]] = []
    bootstrap_rows: list[dict[str, Any]] = []
    staged_models: dict[str, dict[str, Any]] = {}
    outcomes: dict[str, dict[str, Any]] = {}

    for category in config.categories:
        train_positive = records[category]["train"]
        validation_positive = records[category]["validation"]
        test_positive = records[category]["test"]
        train_background = sample_background_by_category(
            background_by_category(records, category, "train"),
            config.max_all_other_negatives,
            f"{seed}|{category}|train_background",
        )
        validation_background = background_by_category(records, category, "validation")
        test_background = background_by_category(records, category, "test")
        train_negative = _sample_records(
            flatten_background(train_background),
            config.max_all_other_negatives,
            f"{seed}|{category}|train",
        )
        validation_regimes = build_negative_regimes(
            validation_background,
            train_positive,
            config,
            f"{seed}|{category}|validation",
        )
        test_regimes = build_negative_regimes(
            test_background,
            train_positive,
            config,
            f"{seed}|{category}|test",
        )
        if any(len(rows) < config.min_split_size for rows in validation_regimes.values()):
            raise ValueError(f"Validation negative regime too small for {category}, seed {seed}.")
        if any(len(rows) < config.min_split_size for rows in test_regimes.values()):
            raise ValueError(f"Test negative regime too small for {category}, seed {seed}.")

        selected_props, ranges, property_statistics = fit_property_component(
            train_positive,
            train_negative,
            config,
        )
        validation_category_rows: dict[str, Sequence[MoleculeRecord]] = {
            category: validation_positive,
            **validation_background,
        }
        test_category_rows: dict[str, Sequence[MoleculeRecord]] = {
            category: test_positive,
            **test_background,
        }

        models: dict[str, FrozenModel] = {}
        selected_by_method: dict[str, tuple[dict[str, str], dict[str, float]]] = {}
        for method in STRUCTURAL_METHODS:
            specs = discover_candidate_specs(train_positive, method, config)
            assessments, patterns, weights = assess_and_select_patterns(
                specs,
                train_positive,
                train_background,
                config,
            )
            selected_by_method[method] = (patterns, weights)
            for assessment in assessments:
                pattern_rows.append(
                    {
                        "seed": seed,
                        "category": category,
                        "method": method,
                        "pattern_id": assessment.spec.pattern_id,
                        "origin": assessment.spec.origin,
                        "source_unit": assessment.spec.source_unit,
                        "smarts": assessment.spec.smarts,
                        "discovery_split": "train",
                        "discovery_count": assessment.spec.discovery_count,
                        "positive_count": assessment.positive_count,
                        "positive_rate": assessment.positive_rate,
                        "negative_rate": assessment.negative_rate,
                        "max_other_rate": assessment.max_other_rate,
                        "category_breadth": assessment.category_breadth,
                        "enrichment": assessment.enrichment,
                        "specificity_gap": assessment.specificity_gap,
                        "quality": assessment.quality,
                        "selected": int(assessment.selected),
                        "rejection_reason": assessment.rejection_reason,
                        "weight": assessment.weight,
                    }
                )
            models[method] = tune_frozen_model(
                category,
                method,
                selected_props,
                ranges,
                patterns,
                weights,
                validation_positive,
                validation_regimes,
                validation_category_rows,
                config.weight_grid,
            )

        models["property_only"] = tune_frozen_model(
            category,
            "property_only",
            selected_props,
            ranges,
            {},
            {},
            validation_positive,
            validation_regimes,
            validation_category_rows,
            (1.0,),
        )
        hybrid_patterns, hybrid_weights = selected_by_method["hybrid"]
        models["hybrid_structure_only"] = tune_frozen_model(
            category,
            "hybrid_structure_only",
            selected_props,
            ranges,
            hybrid_patterns,
            hybrid_weights,
            validation_positive,
            validation_regimes,
            validation_category_rows,
            (0.0,),
        )

        candidate_method = max(
            CANDIDATE_METHODS,
            key=lambda name: (models[name].validation_selection_score, name),
        )
        overall_method = max(
            models,
            key=lambda name: (models[name].validation_selection_score, name),
        )
        test_metrics: dict[str, dict[str, dict[str, float]]] = {}
        test_scores: dict[str, dict[str, np.ndarray]] = {}
        test_cross: dict[str, dict[str, dict[str, float]]] = {}
        for method, model in sorted(models.items()):
            method_metrics, method_scores = evaluate_frozen_model(model, test_positive, test_regimes)
            test_metrics[method] = method_metrics
            test_scores[method] = method_scores
            test_cross[method] = cross_category_rates(model, test_category_rows)
            for regime, metrics in method_metrics.items():
                metric_rows.append(
                    {
                        "seed": seed,
                        "category": category,
                        "method": method,
                        "regime": regime,
                        "split": "test",
                        "threshold_source": "validation_all_other",
                        "threshold": model.threshold,
                        "best_w": model.best_w,
                        "pattern_count": len(model.patterns),
                        "property_count": len(model.selected_props),
                        "validation_selection_score": model.validation_selection_score,
                        **metrics,
                    }
                )
            for evaluated_category, rates in test_cross[method].items():
                cross_rows.append(
                    {
                        "seed": seed,
                        "model_category": category,
                        "method": method,
                        "evaluated_category": evaluated_category,
                        "split": "test",
                        "overlap_policy": "target_exact_members_excluded",
                        **rates,
                    }
                )

        candidate_model = models[candidate_method]
        staged_models[category] = _model_to_json(
            candidate_model,
            seed,
            test_metrics[candidate_method],
            split_digest,
        )

        for baseline_method in ("fixed_smarts", "property_only"):
            for regime in EVALUATION_REGIMES:
                intervals = paired_bootstrap_deltas(
                    test_scores[candidate_method],
                    candidate_model.threshold,
                    test_scores[baseline_method],
                    models[baseline_method].threshold,
                    regime,
                    config.bootstrap_replicates,
                    int(stable_hash(f"{seed}|{category}|{baseline_method}|{regime}")[:8], 16),
                )
                candidate_metrics = test_metrics[candidate_method][regime]
                baseline_metrics = test_metrics[baseline_method][regime]
                bootstrap_rows.append(
                    {
                        "seed": seed,
                        "category": category,
                        "candidate_method": candidate_method,
                        "baseline_method": baseline_method,
                        "regime": regime,
                        "delta_auc": candidate_metrics["auc"] - baseline_metrics["auc"],
                        "delta_balanced_accuracy": candidate_metrics["balanced_accuracy"]
                        - baseline_metrics["balanced_accuracy"],
                        **intervals,
                    }
                )

        candidate_cross = test_cross[candidate_method]
        target_rate = candidate_cross[category]["positive_rate"]
        off_target_rates = [
            rates["positive_rate"]
            for other_category, rates in candidate_cross.items()
            if other_category != category
        ]
        outcomes[category] = {
            "seed": seed,
            "category": category,
            "candidate_method": candidate_method,
            "overall_validation_winner": overall_method,
            "candidate_won_validation": candidate_method == overall_method,
            "pattern_smarts": sorted(candidate_model.patterns.values()),
            "pattern_count": len(candidate_model.patterns),
            "best_w": candidate_model.best_w,
            "test_metrics": test_metrics[candidate_method],
            "baseline_metrics": {
                baseline: test_metrics[baseline] for baseline in ("fixed_smarts", "property_only")
            },
            "target_positive_rate": target_rate,
            "max_off_target_rate": max(off_target_rates, default=0.0),
            "specificity_gap": target_rate - max(off_target_rates, default=0.0),
            "selected_props": list(selected_props),
            "property_ks": property_statistics,
        }

    return {
        "seed": seed,
        "split_rows": split_rows,
        "split_digest": split_digest,
        "metric_rows": metric_rows,
        "pattern_rows": pattern_rows,
        "cross_rows": cross_rows,
        "bootstrap_rows": bootstrap_rows,
        "staged_models": staged_models,
        "outcomes": outcomes,
    }


def _mean(values: Sequence[float]) -> float:
    return float(sum(values) / len(values)) if values else 0.0


def _pattern_stability(pattern_sets: Sequence[set[str]]) -> float:
    if len(pattern_sets) < 2:
        return 1.0
    values: list[float] = []
    for left_index, left in enumerate(pattern_sets):
        for right in pattern_sets[left_index + 1 :]:
            union = left | right
            values.append(len(left & right) / len(union) if union else 1.0)
    return _mean(values)


def build_promotion_rows(
    seed_results: Sequence[dict[str, Any]],
    config: ValidationConfig,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    all_bootstrap = [row for result in seed_results for row in result["bootstrap_rows"]]
    for category in config.categories:
        outcomes = [result["outcomes"][category] for result in seed_results]
        method_counts = Counter(outcome["candidate_method"] for outcome in outcomes)
        dominant_method, dominant_count = sorted(
            method_counts.items(),
            key=lambda item: (-item[1], item[0]),
        )[0]
        method_consistency = dominant_count / len(outcomes)
        pattern_stability = _pattern_stability(
            [set(outcome["pattern_smarts"]) for outcome in outcomes]
        )
        candidate_metrics = [
            metrics
            for outcome in outcomes
            for metrics in outcome["test_metrics"].values()
        ]
        worst_auc = min(metrics["auc"] for metrics in candidate_metrics)
        worst_balanced_accuracy = min(
            metrics["balanced_accuracy"] for metrics in candidate_metrics
        )

        deltas: dict[str, list[float]] = {
            "fixed_auc": [],
            "fixed_ba": [],
            "property_auc": [],
            "property_ba": [],
        }
        for outcome in outcomes:
            for regime in EVALUATION_REGIMES:
                candidate = outcome["test_metrics"][regime]
                fixed = outcome["baseline_metrics"]["fixed_smarts"][regime]
                prop = outcome["baseline_metrics"]["property_only"][regime]
                deltas["fixed_auc"].append(candidate["auc"] - fixed["auc"])
                deltas["fixed_ba"].append(
                    candidate["balanced_accuracy"] - fixed["balanced_accuracy"]
                )
                deltas["property_auc"].append(candidate["auc"] - prop["auc"])
                deltas["property_ba"].append(
                    candidate["balanced_accuracy"] - prop["balanced_accuracy"]
                )

        fixed_bootstrap = [
            row
            for row in all_bootstrap
            if row["category"] == category and row["baseline_method"] == "fixed_smarts"
        ]
        min_auc_ci = min(
            (row["delta_auc_lower"] for row in fixed_bootstrap),
            default=float("nan"),
        )
        min_ba_ci = min(
            (row["delta_balanced_accuracy_lower"] for row in fixed_bootstrap),
            default=float("nan"),
        )
        mean_specificity_gap = _mean(
            [outcome["specificity_gap"] for outcome in outcomes]
        )
        worst_max_off_target = max(
            outcome["max_off_target_rate"] for outcome in outcomes
        )

        failures: list[str] = []
        if len(outcomes) < config.minimum_promotion_seeds:
            failures.append("insufficient_scaffold_seeds")
        if not all(outcome["candidate_won_validation"] for outcome in outcomes):
            failures.append("candidate_not_validation_winner")
        if not all(outcome["pattern_count"] > 0 and outcome["best_w"] < 1.0 for outcome in outcomes):
            failures.append("no_stable_structural_contribution")
        if worst_auc < config.promotion_min_worst_auc:
            failures.append("worst_regime_auc_floor")
        if worst_balanced_accuracy < config.promotion_min_worst_balanced_accuracy:
            failures.append("worst_regime_balanced_accuracy_floor")
        if _mean(deltas["fixed_auc"]) < config.promotion_min_delta_auc:
            failures.append("auc_delta_vs_fixed")
        if _mean(deltas["fixed_ba"]) < config.promotion_min_delta_balanced_accuracy:
            failures.append("balanced_accuracy_delta_vs_fixed")
        if _mean(deltas["property_auc"]) < config.promotion_min_delta_auc:
            failures.append("auc_delta_vs_property_only")
        if _mean(deltas["property_ba"]) < config.promotion_min_delta_balanced_accuracy:
            failures.append("balanced_accuracy_delta_vs_property_only")
        if config.promotion_require_positive_ci:
            if not math.isfinite(min_auc_ci) or min_auc_ci <= 0.0:
                failures.append("paired_auc_ci_not_positive")
            if not math.isfinite(min_ba_ci) or min_ba_ci <= 0.0:
                failures.append("paired_balanced_accuracy_ci_not_positive")
        if mean_specificity_gap < config.promotion_min_specificity_gap:
            failures.append("cross_category_specificity_gap")
        if worst_max_off_target > config.promotion_max_off_target_rate:
            failures.append("cross_category_off_target_rate")
        if method_consistency < config.promotion_min_method_consistency:
            failures.append("candidate_method_instability")
        if pattern_stability < config.promotion_min_pattern_stability:
            failures.append("selected_pattern_instability")

        rows.append(
            {
                "category": category,
                "seed_count": len(outcomes),
                "dominant_candidate_method": dominant_method,
                "method_consistency": method_consistency,
                "pattern_jaccard_stability": pattern_stability,
                "worst_regime_auc": worst_auc,
                "worst_regime_balanced_accuracy": worst_balanced_accuracy,
                "mean_delta_auc_vs_fixed": _mean(deltas["fixed_auc"]),
                "mean_delta_balanced_accuracy_vs_fixed": _mean(deltas["fixed_ba"]),
                "mean_delta_auc_vs_property_only": _mean(deltas["property_auc"]),
                "mean_delta_balanced_accuracy_vs_property_only": _mean(deltas["property_ba"]),
                "minimum_paired_auc_delta_ci_lower_vs_fixed": min_auc_ci,
                "minimum_paired_balanced_accuracy_delta_ci_lower_vs_fixed": min_ba_ci,
                "mean_cross_category_specificity_gap": mean_specificity_gap,
                "worst_max_off_target_positive_rate": worst_max_off_target,
                "promotion_passed": int(not failures),
                "promotion_decision": "eligible_for_manual_review" if not failures else "do_not_promote",
                "failure_reasons": ";".join(failures),
                "automatic_release_performed": 0,
            }
        )
    return rows


def build_pattern_overlap_rows(seed_results: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in seed_results:
        seed = result["seed"]
        outcomes = result["outcomes"]
        categories = sorted(outcomes)
        for left_index, left_category in enumerate(categories):
            left_patterns = set(outcomes[left_category]["pattern_smarts"])
            for right_category in categories[left_index + 1 :]:
                right_patterns = set(outcomes[right_category]["pattern_smarts"])
                union = left_patterns | right_patterns
                rows.append(
                    {
                        "seed": seed,
                        "left_category": left_category,
                        "left_method": outcomes[left_category]["candidate_method"],
                        "right_category": right_category,
                        "right_method": outcomes[right_category]["candidate_method"],
                        "shared_pattern_count": len(left_patterns & right_patterns),
                        "union_pattern_count": len(union),
                        "jaccard": len(left_patterns & right_patterns) / len(union) if union else 1.0,
                    }
                )
    return rows


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _csv_value(value: Any) -> Any:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return format(value, ".12g")
    return value


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(row.get(key, "")) for key in fieldnames})


def _config_payload(config: ValidationConfig) -> dict[str, Any]:
    return {
        field_name: (
            str(value.resolve())
            if isinstance(value, Path)
            else list(value)
            if isinstance(value, tuple)
            else value
        )
        for field_name, value in vars(config).items()
    }


def _implementation_digest() -> str:
    pattern_library_path = APP_DIR / "structural_pattern_library.py"
    digest = hashlib.sha256()
    digest.update(Path(__file__).read_bytes())
    digest.update(pattern_library_path.read_bytes())
    return digest.hexdigest()


def _run_signature(
    config: ValidationConfig,
    input_manifest: dict[str, Any],
    implementation_digest: str,
) -> str:
    config_payload = _config_payload(config)
    config_payload.pop("output_dir", None)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "implementation_sha256": implementation_digest,
        "config": config_payload,
        "input_files": [
            {"category": row["category"], "sha256": row["sha256"]}
            for row in input_manifest["files"]
        ],
    }
    return stable_hash(json.dumps(payload, sort_keys=True, separators=(",", ":")))


def _guard_output_directory(output_dir: Path, run_signature: str) -> None:
    if not output_dir.exists():
        return
    manifest_path = output_dir / "manifest.json"
    existing_files = any(output_dir.iterdir())
    if not existing_files:
        return
    if not manifest_path.exists():
        raise FileExistsError(
            f"Output directory is non-empty but has no validation manifest: {output_dir}. "
            "Choose a new --output-dir."
        )
    existing = json.loads(manifest_path.read_text(encoding="utf-8"))
    if existing.get("run_signature") != run_signature:
        raise FileExistsError(
            f"Output directory belongs to a different validation run: {output_dir}. "
            "Choose a new --output-dir so staged artifacts cannot be mixed."
        )


def run_validation(config: ValidationConfig) -> dict[str, Any]:
    config.validate()
    categories, input_manifest = load_category_inputs(config)
    implementation_digest = _implementation_digest()
    run_signature = _run_signature(config, input_manifest, implementation_digest)
    _guard_output_directory(config.output_dir, run_signature)
    molecule_cache = prepare_molecule_cache(categories)
    seed_results = [run_seed(categories, seed, config, molecule_cache) for seed in config.seeds]
    promotion_rows = build_promotion_rows(seed_results, config)
    overlap_rows = build_pattern_overlap_rows(seed_results)

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    split_rows = [row for result in seed_results for row in result["split_rows"]]
    metric_rows = [row for result in seed_results for row in result["metric_rows"]]
    pattern_rows = [row for result in seed_results for row in result["pattern_rows"]]
    cross_rows = [row for result in seed_results for row in result["cross_rows"]]
    bootstrap_rows = [row for result in seed_results for row in result["bootstrap_rows"]]

    _write_csv(output_dir / "split_assignments.csv", split_rows)
    _write_csv(output_dir / "test_metrics.csv", metric_rows)
    _write_csv(output_dir / "pattern_candidates.csv", pattern_rows)
    _write_csv(output_dir / "cross_category_rates.csv", cross_rows)
    _write_csv(output_dir / "paired_bootstrap_deltas.csv", bootstrap_rows)
    _write_csv(output_dir / "candidate_pattern_overlap.csv", overlap_rows)
    _write_csv(output_dir / "promotion_decisions.csv", promotion_rows)
    for result in seed_results:
        for category, payload in result["staged_models"].items():
            _write_json(output_dir / "models" / f"seed_{result['seed']}" / f"{category}.json", payload)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "run_signature": run_signature,
        "implementation_sha256": implementation_digest,
        "rdkit_version": rdBase.rdkitVersion,
        "config": _config_payload(config),
        "inputs": input_manifest,
        "split_policy": {
            "scope": "global_across_all_product_use_categories",
            "group": "exact_bemis_murcko_with_acyclic_generic_topology_or_canonical_fallback",
            "assignment": "sha256(seed|scaffold_key)",
            "duplicate_policy": "canonical duplicates share a split; target members are excluded from its negatives",
        },
        "leakage_controls": {
            "candidate_discovery": "train_positive_only",
            "property_and_pattern_fit": "train_only",
            "method_weight_threshold_selection": "validation_only",
            "test_threshold": "frozen_validation_all_other_threshold",
            "test_used_for_selection": False,
        },
        "negative_regimes": list(EVALUATION_REGIMES),
        "negative_regime_policy": {
            "training_background": "equal per-category deterministic sample capped by max_all_other_negatives",
            "all_other": "deterministic canonical-structure sample",
            "related_hard": "highest Morgan-Tanimoto similarity to a deterministic training-positive reference set",
            "property_matched": "smallest robust-scaled descriptor distance to the training-positive median",
        },
        "category_semantics": "independent overlapping product-use support; not mutually exclusive classes",
        "excluded_from_comparison": {
            "endocrine_disruptors": "auxiliary hazard endpoint, not a product-use category"
        },
        "deployment": {
            "research_only": True,
            "bundled_models_modified": False,
            "automatic_promotion": False,
            "promotion_pass_means": "eligible_for_manual_review_only",
        },
        "artifacts": [
            "split_assignments.csv",
            "test_metrics.csv",
            "pattern_candidates.csv",
            "cross_category_rates.csv",
            "paired_bootstrap_deltas.csv",
            "candidate_pattern_overlap.csv",
            "promotion_decisions.csv",
            "models/seed_<seed>/<category>.json",
        ],
    }
    _write_json(output_dir / "manifest.json", manifest)
    return {
        "output_dir": str(output_dir.resolve()),
        "seed_count": len(config.seeds),
        "category_count": len(config.categories),
        "promotion_pass_count": sum(row["promotion_passed"] for row in promotion_rows),
        "promotion_rows": promotion_rows,
    }


def _parse_int_tuple(value: str) -> tuple[int, ...]:
    parsed = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    if not parsed:
        raise argparse.ArgumentTypeError("Expected a comma-separated list of integers.")
    return parsed


def _parse_float_tuple(value: str) -> tuple[float, ...]:
    parsed = tuple(float(item.strip()) for item in value.split(",") if item.strip())
    if not parsed:
        raise argparse.ArgumentTypeError("Expected a comma-separated list of numbers.")
    return parsed


def _parse_categories(value: str) -> tuple[str, ...]:
    parsed = tuple(item.strip() for item in value.split(",") if item.strip())
    if not parsed:
        raise argparse.ArgumentTypeError("Expected a comma-separated list of categories.")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Mine fixed/Murcko/BRICS/hybrid structural candidates on training data and "
            "validate frozen category scorers on global scaffold-disjoint holdouts."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--seeds", type=_parse_int_tuple, default=(11, 23, 37))
    parser.add_argument("--categories", type=_parse_categories, default=PRODUCT_USE_CATEGORIES)
    parser.add_argument("--train-fraction", type=float, default=0.60)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument(
        "--limit-per-category",
        type=int,
        default=None,
        help="Deterministic hash sample for bounded pilots; omit for the full local inputs.",
    )
    parser.add_argument("--candidate-pool-size", type=int, default=96)
    parser.add_argument("--max-patterns", type=int, default=24)
    parser.add_argument("--min-fragment-atoms", type=int, default=3)
    parser.add_argument("--max-fragment-atoms", type=int, default=18)
    parser.add_argument("--min-positive-count", type=int, default=3)
    parser.add_argument("--min-positive-prevalence", type=float, default=0.01)
    parser.add_argument("--min-enrichment", type=float, default=1.25)
    parser.add_argument("--min-specificity-gap", type=float, default=0.01)
    parser.add_argument("--min-cross-category-prevalence", type=float, default=0.01)
    parser.add_argument("--max-category-breadth", type=int, default=6)
    parser.add_argument("--property-ks-threshold", type=float, default=0.15)
    parser.add_argument("--weight-grid", type=_parse_float_tuple, default=(0.0, 0.25, 0.5, 0.75, 1.0))
    parser.add_argument("--max-all-other-negatives", type=int, default=2400)
    parser.add_argument("--hard-negative-limit", type=int, default=600)
    parser.add_argument("--bootstrap-replicates", type=int, default=300)
    parser.add_argument("--minimum-promotion-seeds", type=int, default=3)
    parser.add_argument(
        "--allow-nonpositive-bootstrap-ci",
        action="store_true",
        help="Diagnostic only: do not require paired 95%% CI lower bounds above zero.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = ValidationConfig(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        seeds=args.seeds,
        categories=args.categories,
        limit_per_category=args.limit_per_category,
        train_fraction=args.train_fraction,
        validation_fraction=args.validation_fraction,
        candidate_pool_size=args.candidate_pool_size,
        max_patterns=args.max_patterns,
        min_fragment_atoms=args.min_fragment_atoms,
        max_fragment_atoms=args.max_fragment_atoms,
        min_positive_count=args.min_positive_count,
        min_positive_prevalence=args.min_positive_prevalence,
        min_enrichment=args.min_enrichment,
        min_specificity_gap=args.min_specificity_gap,
        min_cross_category_prevalence=args.min_cross_category_prevalence,
        max_category_breadth=args.max_category_breadth,
        property_ks_threshold=args.property_ks_threshold,
        weight_grid=args.weight_grid,
        max_all_other_negatives=args.max_all_other_negatives,
        hard_negative_limit=args.hard_negative_limit,
        bootstrap_replicates=args.bootstrap_replicates,
        minimum_promotion_seeds=args.minimum_promotion_seeds,
        promotion_require_positive_ci=not args.allow_nonpositive_bootstrap_ci,
    )
    summary = run_validation(config)
    print(json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
