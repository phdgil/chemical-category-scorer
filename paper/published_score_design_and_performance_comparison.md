# Published-score design and performance comparison

## Why this note exists
The manuscript now distinguishes between:
- universal QED comparison for every final category, and
- direct published same-category score comparison where a published continuous score actually exists.

In the current project pass, the only confirmed direct published continuous same-category comparator beyond QED is the pesticide-likeness framework of Avram et al. 2014.

## Category with direct published score comparator
- **Pesticides**
  - Published score family: **QEH, QEI, QEF, QEPmax, QEPavg**
  - Source: Avram S et al. *Quantitative estimation of pesticide-likeness for agrochemical discovery*. J Cheminform. 2014;6:42. doi:10.1186/s13321-014-0042-6.

## Design comparison

### Published QEP design (Avram 2014)
- Data-driven? **Yes**, but in a descriptor-distribution fitting sense.
- Positive source: marketed herbicides, insecticides, and fungicides assembled from pesticide references.
- Core descriptors: **MW, LogP, HBA, HBD, RB, aromatic rings**.
- Model form:
  - descriptor-specific desirability functions fitted separately for herbicides, insecticides, fungicides
  - class scores combined by geometric mean
  - overall pesticide score fused as **QEPmax** or **QEPavg**
- Structural patterns: **No explicit structural-pattern component**.
- Negative design in the publication:
  - evaluated against AgroSAR patent pesticides and large random PubChem decoy sets
  - not built around the current retained cross-category PubChem negative doctrine.

### Present study design
- Data-driven? **Yes**, in a broader whole-category classification sense.
- Positive source: full PubChem pesticide-category positives after CID-to-SMILES conversion and deduplication.
- Negative source: other final PubChem categories, with cross-category overlap retained except for exact target overlap, duplicate negatives, and near-positive Tanimoto matches.
- Core descriptors: **MW, FCsp3, TPSA, HBA, HBD, RotBonds, AromaticRings**.
- Structural patterns: **Yes**; final pesticide scorer uses **chloroaromatic, CF3, triazole, nitro, carbamate, urea**.
- Optimization: descriptor-range selection + pattern enrichment + best mixing weight selection.

## Performance comparison

### Published paper's own reported AUCs (Avram 2014, Supplementary Table S7)
- QEH: **0.721 ± 0.007**
- QEI: **0.668 ± 0.003**
- QEF: **0.677 ± 0.003**
- QEPmax: **0.643 ± 0.002**
- QEPavg: **0.651 ± 0.002**

These are the publication's original AgroSAR/decoy-benchmark values and are not directly interchangeable with the current PubChem retained-negative regime.

### Same-regime reimplementation on the present PubChem pesticide task
Source file: `results/final_category_rebuild/pesticide_published_score_comparison.csv`

| Score | AUC | Balanced accuracy |
| --- | ---: | ---: |
| QEH | 0.6644 | 0.6187 |
| QEI | 0.6587 | 0.6520 |
| QEF | 0.6694 | 0.6704 |
| QEPmax | 0.6811 | 0.6559 |
| QEPavg | 0.6734 | 0.6429 |
| Final rebuilt pesticide scorer (this study) | 0.7527 | 0.6840 |

## Interpretation
- The published pesticide score is also data-driven, so the novelty claim for pesticides must not be "first data-driven pesticide score".
- The actual differentiator of the present study is the **broader category-reconstruction doctrine**:
  - full positive PubChem rebuild,
  - retained cross-category hard negatives,
  - explicit structural-pattern evidence,
  - same-regime benchmark against both **QED** and the published **QEP** family.
- On the present regime, the final rebuilt pesticide scorer outperforms the best published fusion comparator (**QEPmax**) by:
  - **+0.0716 AUC**
  - **+0.0281 balanced accuracy**

## Important implementation note
- The present same-regime QEP comparison is a transparent reimplementation from the published supplementary coefficients.
- Avram 2014 used ChemAxon descriptor generation; the current project uses RDKit descriptor calculation.
- Therefore this is a **same-regime comparator reimplementation**, not a byte-identical replay of the historical software environment.
