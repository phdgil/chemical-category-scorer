# Chemical Category Scorer

Local-first chemical product category scoring from molecular structure.

This repository ships the project in two delivery modes:
1. **Desktop app** for single-molecule and batch CSV scoring on a local workstation
2. **Python library** for direct import, similar to using an RDKit scoring helper

The local-first design matters for industrial chemistry use because molecular structures do not need to be uploaded to an external service.

## Article-synchronized release

Version **2.1.0** is synchronized with the scoring functions reported in the associated manuscript. The desktop app and Python library load the same four bundled model definitions, thresholds, and cross-category calibration:

- `han_endocrine_disruptors`
- `final_flavor_fragrance`
- `final_pesticides`
- `final_surfactants`

Flavoring agents and fragrances are represented by one merged score because their original positive sets shared 1,071 exact structures and the merged function passed three held-out evaluations. Candidate functions for animal drugs, human drugs, cosmetics, food additives, food-contact substances, fragrances and flavoring agents separately, and solvents are not part of the public 2.1.0 panel. The pesticide function in 2.1.0 adds fifteen atom-neighborhood fragments discovered by positive- versus negative-seeded personalized PageRank on a molecule-fragment network; all fifteen were selected in every held-out fold.

## Desktop app

<img src="docs/desktop_app_opened.png" alt="Chemical Category Scorer desktop app" width="900" />

The desktop entry point is installable from the Python package. After `pip install .`, this command should open the app:

```bash
chemical-category-scorer-desktop
```

The installed command also supports the same model-list smoke check as the source-tree launcher:

```bash
chemical-category-scorer-desktop --list-models
```

### Run from the repository

From the `app` directory:

```bat
run_desktop_app.bat
```

or:

```bash
python app/desktop_app.py
```

### Desktop app capabilities

- single-SMILES scoring
- batch CSV scoring
- local JSON-backed score models
- output without sending structures to a hosted API

The desktop interface accepts one SMILES string or a batch CSV. Select **All article models** to display the complete synchronized score panel. Product-use results and the endocrine-disruption auxiliary signal are reported separately.

## Python library

### Install from the public GitHub repository

```bash
pip install "chemical-category-scorer @ git+https://github.com/phdgil/chemical-category-scorer.git"
```

### Install from a local clone

```bash
pip install .
```

### Python usage

```python
from rdkit import Chem
from chemical_category_scorer import flavor_fragrance, pesticides, details_mol, available_models

mol = Chem.MolFromSmiles("ClC1=C(Cl)C(C#N)=C(Cl)C(C#N)=C1Cl")
score = pesticides(mol)
sensory_score = flavor_fragrance(mol)
info = details_mol(mol, model_id="final_pesticides")
models = available_models()

print(score)
print(sensory_score)
print(info.score, info.threshold, info.decision, info.matched_patterns)
print(models)
```

Use `details_smiles(smiles, model_id=...)` when starting from a SMILES string. It returns the score, threshold, decision, descriptor and structural components, and matched structural patterns.

### Command-line entry points after installation

```bash
chemical-category-scorer --list-models
chemical-category-scorer --score "CCO" --model-id final_pesticides
chemical-category-scorer-desktop
chemical-category-scorer-desktop --list-models
```

## Available scoring categories

Version 2.1.0 exposes exactly the four scoring functions reported in the associated article. The product-use scores and Han endocrine-disruption score answer different questions.

Product-use scorers:

- flavor_fragrance (`final_flavor_fragrance`)
- pesticides
- surfactants

Auxiliary hazard/activity signal:

- endocrine_disruptors (`han_endocrine_disruptors`)

## All-model interpretation

All-model scoring is a multi-label screening view. Each model has its own threshold, fitted independently against its own positive set and constructed background. Scores and margins are model-specific heuristics; they are not calibrated probabilities and are not validated cross-model distances.

Multiple threshold-positive product-use categories can be chemically reasonable because the product-use categories overlap. A representative product-use result is emitted only when exactly one score reaches its high-specificity cross-category threshold; otherwise it remains unresolved. The endocrine-disruption signal is reported separately as an auxiliary Han Se-eum hazard/activity signal and does not replace a product-use category.

## Public-model audit

The reproducible audit entry point is `app.public_model_audit`. Run these commands from the repository root:

```bash
python -m app.public_model_audit audit --scores-out docs/cross_category_probe_scores.csv --summary-out docs/cross_category_probe_summary.csv
python -m app.public_model_audit pattern-overlap --csv-out docs/public_model_pattern_overlap.csv --markdown-out docs/cross_category_pattern_audit.md
python -m app.public_model_audit pattern-candidates --input-csv path/to/category_smiles.csv --category-column category --smiles-column SMILES --output-csv results/pattern_candidates/fixed_murcko_brics_hybrid.csv
```

The candidate input CSV must contain one category column and one SMILES column. The local-only `--use-final-rebuild-inputs` shortcut is for development datasets and is not required by the released package. The seven-molecule probe audit is intended for post-deployment diagnostics on caffeine, aspirin, DDT, bisphenol A, vanillin, SDS, and ethanol. It is not the formal publication benchmark. BRICS, Murcko scaffolds, fixed SMARTS, and hybrid candidates should be compared as a preregistered held-out ablation before any model replacement claim.

## Structural-pattern validation loop

When the local development positive CSVs are available, run the leakage-resistant comparison with:

```bash
python -m app.structural_pattern_validation --seeds 11,23,37
```

For a quick end-to-end check, add `--limit-per-category 80 --bootstrap-replicates 10`. The loop assigns duplicate structures and scaffold groups globally to train, validation, or test; mines fixed SMARTS, Murcko, BRICS, and hybrid candidates from training positives only; freezes the method, mixture weight, and threshold on validation; and reports three untouched test-negative regimes plus a cross-category positive-rate matrix. Results and runtime-compatible candidate JSONs are written under `results/structural_pattern_validation/`.

This command never replaces `app/data/models`. Even a passing gate means only `eligible_for_manual_review`; a candidate that lacks multi-seed stability, held-out improvement over fixed SMARTS and property-only baselines, or cross-category specificity is marked `do_not_promote`.

Use a new `--output-dir` when changing inputs, seeds, or loop settings. The command refuses to mix artifacts from runs with different signatures.

## Repository contents

- `app/` desktop application and scoring engine
- `chemical_category_scorer/` importable Python API
- `docs/` app screenshot and pattern-testing reference files
- `paper/` and `results/` are local ignored working areas unless explicitly published

## Privacy

- no web upload is required for scoring
- inputs are read from local files only
- outputs are written to local files only
- model JSON files are bundled with the repository
