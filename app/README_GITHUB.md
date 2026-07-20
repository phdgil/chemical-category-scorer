# Chemical Category Scorer Desktop App

Local-first desktop application for broad chemical product-category scoring from molecular structure.

## Why a desktop app
Chemical manufacturers often treat molecular structures as confidential business assets. This app is designed so that scoring can run on a local workstation without sending SMILES strings, descriptors, or batch files to an external server. That makes it suitable for GitHub distribution while still preserving offline use inside companies.

## Current capabilities
- single-SMILES scoring
- batch CSV scoring
- local JSON-backed score models
- simple single-molecule and batch scoring workflow focused on the final manuscript models
- advanced scorer rebuilding is available through project scripts, not through the desktop UI
- local-first execution for confidential industrial structures

## Final manuscript category set
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

## Install dependencies
```bash
pip install -r requirements.txt
```

## Verify installation
```bash
python desktop_app.py --self-test
```

## Recommended GitHub release contents
- `desktop_app.py`
- `algorithm_score_engine.py`
- `data/models/`
- `data/evidence_panels/`
- `requirements.txt`
- `run_desktop_app.bat`
- example input/output CSV files

## Repository note
If this app is published on GitHub for public reuse, keep the scoring local-first. Do not replace local scoring with a remote API by default, because the local privacy guarantee is part of the practical value of this project.


## Python library use
After `pip install .` from the repository root, users can import the scorers directly:

```python
from rdkit import Chem
from chemical_category_scorer import pesticides, details_mol

mol = Chem.MolFromSmiles("ClC1=C(Cl)C(C#N)=C(Cl)C(C#N)=C1Cl")
score = pesticides(mol)
info = details_mol(mol, model_id="final_pesticides")
```

This gives the repository two local-first delivery modes: a desktop app and an importable Python library.


Built-in final categories now include endocrine disruptors in both the desktop app and the Python library. The deployed endocrine scorer is the Han Se-eum model selected after cross-regime validation against the uncapped GJC endocrine rebuild on the final broad-category evaluation set.
