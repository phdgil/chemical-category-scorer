# Journal of Cheminformatics manuscript draft template

## Title page
**Proposed title**  
Structure-based chemical product-category scoring from molecular structure: full-positive PubChem broad-category rebuild, QED comparison, and desktop deployment

**Authors**  
[Fill full author list]

**Affiliations**  
[Fill affiliations]

**Corresponding author**  
[Name, email, address]

---

## Abstract
Chemical product-category suggestion from molecular structure is useful for library triage, chemical inventory analysis, and downstream safety or product-screening workflows, but most descriptor-based scoring studies have focused on drug-likeness rather than broader chemical product classes. Here we report a structure-based scoring framework for ten broad chemical product categories derived from the PubChem Classification Browser HID 72 tree: animal drugs, human drugs, cosmetics, flavoring agents, food additives, food contact substances, fragrances, pesticides, solvents, and surfactants. The final run used all available category positives in this project after CID-to-SMILES conversion and SMILES deduplication, retained cross-category negatives, and combined positive-set physicochemical property ranges with structural patterns in one interpretable scorer. All final scorers were benchmarked against raw RDKit quantitative estimate of drug-likeness (QED) on the same rebuilt evaluation sets. The strongest categories were surfactants (AUC 0.9778), flavoring agents (0.9389), solvents (0.9295), and fragrances (0.9254), while harder broad categories such as food contact substances (0.6213), cosmetics (0.6179), and animal drugs (0.6136) remained above QED. The final framework was also translated into a desktop application for practical scoring and batch use. These results show that broad chemical product-category scoring beyond classic drug-likeness is technically feasible and that same-regime comparison against QED provides a clear minimum benchmark for publication and deployment.

### Scientific Contribution
We extend QED-style interpretable molecular scoring from drug-likeness toward ten broad PubChem chemical product categories using full-positive uncapped category rebuilds. We show that every final category scorer outperforms raw RDKit QED on the same retained-cross-category evaluation set. We also provide a deployable desktop implementation that turns the scoring framework into a reusable cheminformatics workflow.

**Keywords:** cheminformatics; molecular scoring; QED; PubChem Classification Browser; structural alerts; chemical product categories; descriptor-based scoring; interpretable models

---

## Introduction
[Draft this section around:]
- QED as the best-known descriptor-based molecular desirability score.
- Why chemical product-category suggestion is broader than drug-likeness.
- Why PubChem HID 72 provides a scalable category source.
- Why broad categories such as cosmetics and food contact substances are hard and worth testing explicitly.
- Study objective: build broad-category scorers that beat same-regime QED.

## Methods
### Data source and category definition
- PubChem Classification Browser HID 72 chemical classes.
- Final categories analyzed: animal drugs, human drugs, cosmetics, flavoring agents, food additives, food contact substances, fragrances, pesticides, solvents, surfactants.
- Positive count definition: unique SMILES after CID-to-SMILES conversion and deduplication.

### Positive and negative set construction
- Positives: all available category positives in this run.
- Negatives: pooled from the other final categories.
- Retained cross-category overlap policy.
- Removal only for exact target-positive overlap, duplicate negative SMILES, and near-positive Tanimoto matches.

### Scoring framework
- descriptor ranges from the positive set
- structural-pattern selection
- combined weighted score
- Lee-style Bayesian optimization when it outperforms grid search
- category priors when available

### Baseline comparison
- Raw RDKit QED on the same rebuilt evaluation set.
- Report AUC, balanced accuracy, and delta versus QED.

### Software implementation
- desktop application
- batch CSV scoring
- model packaging

## Results
### Final full-positive category counts
Use `final_category_rebuild_summary.csv` and `full_positive_refresh_counts.csv`.

### Final scorer performance versus QED
Main table source: `final_category_rebuild_qed_comparison.csv`

Suggested narrative:
- surfactants, flavoring agents, solvents, and fragrances are strongest.
- food additives and human drugs are moderate but clearly above QED.
- pesticides improves after removing the positive cap and still beats QED.
- food contact substances, cosmetics, and animal drugs remain harder broad categories but still exceed QED in the final rebuilt regime.

### Structural-pattern content
Report selected property counts and pattern counts from `final_category_rebuild_summary.csv`.

### Regime context
Use `regime_separated_comparison_table.csv` as supplement/context only.

## Discussion
[Discuss:]
- Why full-positive rebuild mattered.
- Why QED is a strong but incomplete baseline.
- Why some broad categories remain difficult despite beating QED.
- Why retained cross-category negatives were scientifically preferable to aggressive overlap removal.
- Why the framework is useful as an interpretable cheminformatics scoring approach.

## Conclusions
- Broad product-category scorers beyond classic drug-likeness are feasible.
- Same-regime QED comparison gives a clean publication gate.
- The full-positive rebuild strengthens the manuscript materially.
- The desktop implementation demonstrates practical utility.

## List of abbreviations
QED, quantitative estimate of drug-likeness  
AUC, area under the receiver operating characteristic curve  
KS, Kolmogorov-Smirnov statistic  
HID, hierarchy identifier  
SMILES, simplified molecular-input line-entry system  
QSAR, quantitative structure-activity relationship

## Declarations
### Availability of data and materials
All project tables, rebuilt models, and manuscript-preparation artifacts are available under the project result and paper folders. Add repository/archive link here before submission.

### Competing interests
The authors declare that they have no competing interests.

### Funding
[Fill]

### Authors’ contributions
[Fill]

### Acknowledgements
[Fill]

### Authors’ information
[Optional]

## References
Use Basic Springer reference style.

Seed references to include:
1. Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL (2012) Quantifying the chemical beauty of drugs. Nat Chem 4:90-98. doi:10.1038/nchem.1243.
2. Ertl P, Schuffenhauer A (2009) Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. J Cheminform 1:8. doi:10.1186/1758-2946-1-8.
3. Sorokina M, Steinbeck C (2019) NaPLeS: a natural products likeness scorer-web application and database. J Cheminform 11:55. doi:10.1186/s13321-019-0378-z.
4. Cortes-Ciriano I (2016) Bioalerts: a python library for the derivation of structural alerts from bioactivity and toxicity data sets. J Cheminform 8:13. doi:10.1186/s13321-016-0125-7.
5. Add PubChem Classification Browser citation/source and RDKit citation.
