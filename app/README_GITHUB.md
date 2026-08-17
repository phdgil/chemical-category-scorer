# Chemical Category Scorer Desktop App

Local-first desktop application for broad chemical product-category scoring from molecular structure.

Version **2.1.0** is synchronized with the associated manuscript and Python library. All three interfaces—the source application, packaged desktop archive, and importable library—use the same four model JSON files, thresholds, and cross-category calibration.

## Why a desktop app
Chemical manufacturers often treat molecular structures as confidential business assets. This app is designed so that scoring can run on a local workstation without sending SMILES strings, descriptors, or batch files to an external server. That makes it suitable for GitHub distribution while still preserving offline use inside companies.

## Current capabilities
- single-SMILES scoring
- batch CSV scoring
- local JSON-backed score models
- simple single-molecule and batch scoring workflow for the deployed category scorers
- all-public-model scoring for transparent multi-label screening
- local-first execution for confidential industrial structures

## Available categories

Version 2.1.0 contains exactly the four scoring functions reported in the associated article: three product-use functions and one auxiliary Han Se-eum endocrine-disruption hazard/activity function. The pesticide function includes fifteen consensus atom-neighborhood fragments identified by held-out bipartite-network random-walk analysis.

Product-use scorers:

- flavor_fragrance (`final_flavor_fragrance`)
- pesticides
- surfactants

Auxiliary signal:

- endocrine_disruptors (`han_endocrine_disruptors`)

Flavoring agents and fragrances are represented by the merged `flavor_fragrance` score because the original positive sets shared 1,071 exact structures and the merged function passed three held-out evaluations. Categories rejected during exact-overlap-controlled publication screening are not exposed by the release.

## Score interpretation

All-model output should be read as independently thresholded multi-label screening. A molecule can cross more than one product-use threshold because the categories are broad and chemically overlapping. Raw scores and margins are not calibrated probabilities and are not comparable distances across models.

The app now distinguishes **category-enriched evidence** from **shared/nonspecific evidence** using a second threshold calibrated against hard cross-category structures. It reports one representative product-use category only when exactly one product-use score reaches this higher-specificity operating point. Otherwise, the representative result is explicitly unresolved and the complete score panel remains visible.

The endocrine-disruption result is an auxiliary Han signal, not a peer product-use category. Batch and audit reporting should keep this signal separate from the representative product-use category.

## Privacy model
- no web upload is required for scoring
- inputs are read from local files only
- outputs are written to local files only
- model JSON files are bundled with the repository
- users can disconnect from the network after installation and still run the app

## Run the desktop app
From the `app` directory:

```bat
run_desktop_app.bat
```

or:

```bash
python desktop_app.py
```

Enter one SMILES string for interactive scoring, or select an input CSV and SMILES column for batch scoring. Selecting the all-model option displays all three product-use scores together and keeps the endocrine-disruption result in its separate auxiliary-hazard section.

## Install dependencies
```bash
pip install -r requirements.txt
```

From the repository root, the package install exposes both command-line entry points:

```bash
pip install .
chemical-category-scorer --list-models
chemical-category-scorer-desktop --list-models
```

## Verify installation
```bash
python desktop_app.py --self-test
```

From the repository root, these smoke checks verify the package import path and direct launcher path:

```bash
python -c "import app.desktop_app"
python app/desktop_app.py --list-models
```

## Public audit commands

The audit CLI module is `app.public_model_audit`. Run it from the repository root in this shape:

```bash
python -m app.public_model_audit audit --scores-out docs/cross_category_probe_scores.csv --summary-out docs/cross_category_probe_summary.csv
python -m app.public_model_audit pattern-overlap --csv-out docs/public_model_pattern_overlap.csv --markdown-out docs/cross_category_pattern_audit.md
python -m app.public_model_audit pattern-candidates --input-csv path/to/category_smiles.csv --category-column category --smiles-column SMILES --output-csv results/pattern_candidates/fixed_murcko_brics_hybrid.csv
```

The audit command is for the seven post-deployment probe molecules: caffeine, aspirin, DDT, bisphenol A, vanillin, SDS, and ethanol. These probes are diagnostic examples for all-model interpretation, not the formal publication benchmark.

The pattern-candidate input must provide category and SMILES columns. The local-only `--use-final-rebuild-inputs` shortcut is for development datasets and is not required by the release. The command is experimental: any fixed-SMARTS, Murcko, BRICS, or hybrid claim should come from a held-out ablation, not from replacing deployed models by inspection.

## Recommended GitHub release contents
- `desktop_app.py`
- `algorithm_score_engine.py`
- `data/models/han_endocrine_disruptors.json`
- `data/models/final_flavor_fragrance.json`
- `data/models/final_pesticides.json`
- `data/models/final_surfactants.json`
- `data/evidence_panels/`
- `requirements.txt`
- `run_desktop_app.bat`

## Python library use
After `pip install .` from the repository root, users can import the scorers directly:

```python
from rdkit import Chem
from chemical_category_scorer import pesticides, details_mol

mol = Chem.MolFromSmiles("ClC1=C(Cl)C(C#N)=C(Cl)C(C#N)=C1Cl")
score = pesticides(mol)
info = details_mol(mol, model_id="final_pesticides")
```

This repository provides two local-first delivery modes: a desktop app and an importable Python library.
