# Cross-category pattern audit

This report is derived only from the public model JSON configs listed in `app/data/app_release_config.json`; it does not mutate thresholds or model files.

The 65-row `docs/app_test_smiles_by_category.csv` file is a structural pattern unit-test panel. It verifies that configured SMARTS can match fixed probe molecules, but it is not a labeled classification benchmark and should not be used for performance claims.

Public pattern units inspected: 22. Exact canonical SMARTS reused across more than one model: 0.

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

- No exact canonical SMARTS reuse was found.

## Interpretation

Shared fixed SMARTS can help explain why simple molecules can trigger multiple independent category thresholds. The seven-probe audit is exploratory cross-category behavior checking, not held-out validation. A 30-per-category exploratory pilot is only justified when a concrete sampled input panel is recorded; this repository slice does not add such evidence and makes no accuracy or BRICS-improvement claim.
