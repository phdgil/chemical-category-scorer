# Chemical Category Scorer

Local-first chemical product category scoring from molecular structure.

This repository ships the project in two delivery modes:
1. **Desktop app** for single-molecule and batch CSV scoring on a local workstation
2. **Python library** for direct import, similar to using an RDKit scoring helper

The local-first design matters for industrial chemistry use because molecular structures do not need to be uploaded to an external service.

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
from chemical_category_scorer import pesticides, details_mol, available_models

mol = Chem.MolFromSmiles("ClC1=C(Cl)C(C#N)=C(Cl)C(C#N)=C1Cl")
score = pesticides(mol)
info = details_mol(mol, model_id="final_pesticides")
models = available_models()

print(score)
print(info.score, info.threshold, info.decision, info.matched_patterns)
print(models)
```

### Command-line entry points after installation

```bash
chemical-category-scorer --list-models
chemical-category-scorer --score "CCO" --model-id final_pesticides
chemical-category-scorer-desktop
chemical-category-scorer-desktop --list-models
```

## Available scoring categories

The public app exposes ten product-use scorers and one auxiliary hazard/activity signal. Product-use scores and the Han endocrine-disruption score are reported together for convenience, but they answer different questions.

Product-use scorers:

- animal_drugs
- human_drugs
- cosmetics
- flavoring_agents
- food_additives
- food_contact_substances
- fragrances
- pesticides
- solvents
- surfactants

Auxiliary hazard/activity signal:

- endocrine_disruptors (`han_endocrine_disruptors`)

## All-model interpretation

All-model scoring is a multi-label screening view. Each model has its own threshold, fitted independently against its own positive set and constructed background. Scores and margins are model-specific heuristics; they are not calibrated probabilities and are not validated cross-model distances.

Multiple threshold-positive product-use categories can be chemically reasonable because the product-use categories overlap. The highest raw product-use score is a ranked screening suggestion, not proof of the true or closest category. The endocrine-disruption signal is reported separately as an auxiliary Han Se-eum hazard/activity signal and should not silently replace a product-use category.

## Public-model audit

The reproducible audit entry point is `app.public_model_audit`. Run these commands from the repository root:

```bash
python -m app.public_model_audit audit --scores-out docs/cross_category_probe_scores.csv --summary-out docs/cross_category_probe_summary.csv
python -m app.public_model_audit pattern-overlap --csv-out docs/public_model_pattern_overlap.csv --markdown-out docs/cross_category_pattern_audit.md
python -m app.public_model_audit pattern-candidates --input-csv path/to/category_smiles.csv --category-column category --smiles-column SMILES --output-csv results/pattern_candidates/fixed_murcko_brics_hybrid.csv
```

The candidate input CSV must contain one category column and one SMILES column. The local-only `--use-final-rebuild-inputs` shortcut is available when all ten files under `app/output/final_category_rebuild/inputs/` exist; it fails clearly when they do not. The seven-molecule probe audit is intended for post-deployment diagnostics on caffeine, aspirin, DDT, bisphenol A, vanillin, SDS, and ethanol. It is not the formal ten-category classification benchmark. BRICS, Murcko scaffolds, fixed SMARTS, and hybrid candidates should be compared as a preregistered held-out ablation before any model replacement claim.

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
