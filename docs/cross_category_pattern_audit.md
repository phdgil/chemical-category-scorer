# Cross-category pattern audit

This report is derived only from the public model JSON configs listed in `app/data/app_release_config.json`; it does not mutate thresholds or model files.

The 65-row `docs/app_test_smiles_by_category.csv` file is a structural pattern unit-test panel. It verifies that configured SMARTS can match fixed probe molecules, but it is not a labeled classification benchmark and should not be used for performance claims.

Public pattern units inspected: 31. Exact canonical SMARTS reused across more than one model: 14.

## Pattern-unit panel data quality

- Rows: 65
- Columns: 12
- Missing cells: 0
- Exact duplicate rows: 0
- Unique test SMILES: 25
- Invalid test SMILES: 0
- Invalid pattern SMARTS: 0

Repeated probe SMILES are expected in this file because one molecule can be used to test several pattern units across different public models. These repeats are therefore not classification-label duplicates.

## Shared fixed SMARTS

- `C(=O)OC` appears in 7 models: final_animal_drugs; final_cosmetics; final_flavoring_agents; final_food_additives; final_food_contact_substances; final_fragrances; final_human_drugs.
- `CCCCCC` appears in 7 models: final_animal_drugs; final_cosmetics; final_food_contact_substances; final_fragrances; final_human_drugs; final_solvents; final_surfactants.
- `NC(=O)N` appears in 3 models: final_animal_drugs; final_human_drugs; final_pesticides.
- `OC(=O)/C=C/c1ccccc1` appears in 4 models: final_cosmetics; final_flavoring_agents; final_food_additives; final_fragrances.
- `OC(=O)N` appears in 2 models: final_human_drugs; final_pesticides.
- `OCCO` appears in 5 models: final_animal_drugs; final_cosmetics; final_food_contact_substances; final_human_drugs; final_surfactants.
- `Oc1ccccc1` appears in 2 models: final_animal_drugs; final_human_drugs.
- `S(=O)=O` appears in 3 models: final_cosmetics; final_food_contact_substances; final_surfactants.
- `[C&H1]=O` appears in 3 models: final_flavoring_agents; final_food_additives; final_fragrances.
- `[Cl,Br,F,I]` appears in 2 models: final_animal_drugs; final_human_drugs.
- `[N&+]` appears in 3 models: final_cosmetics; final_human_drugs; final_surfactants.
- `[N&+](=O)[O&-]` appears in 2 models: final_human_drugs; final_pesticides.
- `c1ccccc1C(=O)c1ccccc1` appears in 3 models: final_cosmetics; final_food_contact_substances; final_fragrances.
- `c1cnc[n&H1]1` appears in 2 models: final_animal_drugs; final_human_drugs.

## Interpretation

Shared fixed SMARTS can help explain why simple molecules can trigger multiple independent category thresholds. The seven-probe audit is exploratory cross-category behavior checking, not held-out validation. A 30-per-category exploratory pilot is only justified when a concrete sampled input panel is recorded; this repository slice does not add such evidence and makes no accuracy or BRICS-improvement claim.
