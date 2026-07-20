# Category-specific scoring literature map

## Why this file exists
The current manuscript is strongest when it shows exactly where prior category-specific scoring exists, where only adjacent functional-use/classification work exists, and where the present study is filling a real gap.

## Core primary-source references to add

### 1. QED template / closest conceptual ancestor
1. Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. Quantifying the chemical beauty of drugs. *Nature Chemistry*. 2012;4:90-98. doi:10.1038/nchem.1243.
   - Use in Introduction/Discussion as the direct template for descriptor-based continuous scoring.

### 2. Pesticide-specific scoring literature
2. Tice CM. Selecting the right compounds for screening: does Lipinski's Rule of 5 for pharmaceuticals apply to agrochemicals? *Pest Management Science*. 2001;57:3-16. doi:10.1002/1526-4998(200101)57:1<3::AID-PS269>3.0.CO;2-6.
   - Early rule-based pesticide/agrochemical profiling.
3. Hao G, Dong Q, Yang G. A comparative study on the constitutive properties of marketed pesticides. *Molecular Informatics*. 2011;30:614-622. doi:10.1002/minf.201100020.
   - Descriptor-level constitutive-property study of marketed pesticides.
4. Avram S, Funar-Timofei S, Borota A, Chennamaneni SR, Manchala AK, Muresan S. Quantitative estimation of pesticide-likeness for agrochemical discovery. *Journal of Cheminformatics*. 2014;6:42. doi:10.1186/s13321-014-0042-6.
   - Strongest direct precedent for a category-specific continuous score.
   - Introduces QEH/QEI/QEF/QEP from six interpretable descriptors and validates against large agrochemical datasets.

### 3. Fragrance-specific scoring / fragrance-like space
5. Ruddigkeit L, Awale M, Reymond J-L. Expanding the fragrance chemical space for virtual screening. *Journal of Cheminformatics*. 2014;6:27. doi:10.1186/1758-2946-6-27.
   - Defines an explicit fragrance-like property range and validates recovery of fragrance families.
   - Not a continuous fragrance score in the QED sense, but a strong direct precedent for category-focused fragrance-likeness design.

### 4. Food additive / flavoring agents / GRAS literature
6. Sprous DG, Salemme FR. A comparison of the chemical properties of drugs and FEMA/FDA notified GRAS chemical compounds used in the food industry. *Food and Chemical Toxicology*. 2007;45:1419-1427. doi:10.1016/j.fct.2007.02.004.
   - Strong descriptor-space comparison of GRAS compounds versus drugs.
   - Useful for Introduction and Discussion when motivating food-additive/flavoring-agent scoring.
7. Medina-Franco JL, Martínez-Mayorga K, Peppard TL, Del Rio A. Chemoinformatic analysis of GRAS (Generally Recognized as Safe) flavor chemicals and natural products. *PLOS ONE*. 2012;7:e50798. doi:10.1371/journal.pone.0050798.
   - Strong chemoinformatic analysis of FEMA GRAS flavor chemicals.
   - Useful direct precedent for flavoring-agent chemical-space analysis, though not a continuous score.

### 5. Functional-use prediction spanning cosmetics / fragrance / flavorant / surfactant / solvent
8. Isaacs KK, Goldsmith M-R, Egeghy P, Phillips K, Brooks R, Hong T, Wambaugh JF. Characterization and prediction of chemical functions and weight fractions in consumer products. *Toxicology Reports*. 2016;3:723-732. doi:10.1016/j.toxrep.2016.08.011.
   - Uses CosIng-centered personal-care product data, harmonized functions, and QSPR classification models.
   - Important for cosmetics-related functional-use context.
9. Phillips KA, Wambaugh JF, Grulke CM, Dionisio KL, Isaacs KK. High-throughput screening of chemicals as functional substitutes using structure-based classification models. *Green Chemistry*. 2017;19:1063-1074. doi:10.1039/C6GC02744J.
   - Very important adjacent prior art: 41 QSUR functional-use models including fragrance, flavorant, surfactant, and related functions.
   - The paper also shows that broad functional categories such as solvent can be difficult when the category is too heterogeneous.

### 6. Solvent-focused category modeling
10. Gramatica P, Navas N, Todeschini R. Classification of organic solvents and modelling of their physico-chemical properties by chemometric methods using different sets of molecular descriptors. *TrAC Trends in Analytical Chemistry*. 1999;18:461-471. doi:10.1016/S0165-9936(99)00115-6.
   - Strong solvent-category precedent based on molecular descriptors and chemometric classification.
11. Alder CM, Hayler JD, Henderson RK, et al. Toward a more holistic framework for solvent selection. *Organic Process Research & Development*. 2016;20:760-773. doi:10.1021/acs.oprd.6b00015.
   - Not a structure-only score, but a major solvent-classification/selection reference useful for Discussion.

### 7. General category-likeness / fragment-likeness templates
12. Jayaseelan KV, Moreno P, Truszkowski A, Ertl P, Steinbeck C. Natural product-likeness score revisited: an open-source, open-data implementation. *BMC Bioinformatics*. 2012;13:106. doi:10.1186/1471-2105-13-106.
13. Brunner S, Fink L, Reymond J-L. ChEMBL-Likeness Score and Database GDBChEMBL. *Frontiers in Chemistry*. 2020;8:46. doi:10.3389/fchem.2020.00046.
   - These are useful method analogues for fragment-/substructure-frequency scoring logic.
14. Karmaus AL, Filer DL, Martin MT, Houck KA. Evaluation of food-relevant chemicals in the ToxCast high-throughput screening program. *Food and Chemical Toxicology*. 2016;92:188-196. doi:10.1016/j.fct.2016.04.012.
   - Useful adjacent precedent for food contact substances as a chemically intermediate food-relevant category.
15. Huang C, Yang Y, Chen X, Wang C, Li Y, Zheng C, Wang Y. Large-scale cross-species chemogenomic platform proposes a new drug discovery strategy of veterinary drug from herbal medicines. *PLOS ONE*. 2017;12:e0184880. doi:10.1371/journal.pone.0184880.
   - Useful adjacent precedent for animal-drug-oriented drug-likeness based on similarity to approved veterinary drugs.

## Category-by-category prior-art status

| Final manuscript category | Direct continuous score precedent? | Best prior art to cite | Manuscript message |
| --- | --- | --- | --- |
| Pesticides | Yes | Tice 2001; Hao 2011; Avram 2014 | Strongest direct prior-art category. Our work differs by using PubChem broad-category positives, same-regime QED comparison, and descriptor+pattern rebuilding. |
| Fragrances | Partial | Ruddigkeit 2014; Phillips 2017 | There is fragrance-like property-space design and fragrance family recovery, but not the same QED-style broad-category score. |
| Flavoring agents | Partial | Sprous & Salemme 2007; Medina-Franco 2012; Phillips 2017 | Strong chemical-space / functional-use prior art exists, but not the same broad-category continuous score. |
| Food additives | Partial | Sprous & Salemme 2007; Medina-Franco 2012; Isaacs 2016; Phillips 2017 | Adjacent GRAS and functional-use literature exists; direct food-additive continuous scoring appears scarce. |
| Cosmetics | Indirect | Isaacs 2016; Phillips 2017 | Cosmetics has good functional-use / personal-care classification context, but no clear canonical cosmetic-likeness score was found in this pass. |
| Food contact substances | Very limited | Karmaus 2016 plus regulatory/priority literature | Treat as a genuine gap category for continuous scoring, but acknowledge adjacent food-relevant chemical-space clustering work. |
| Solvents | Partial | Gramatica 1999; Alder 2016; Phillips 2017 | Solvent classification exists, but broad solvent use is heterogeneous; this supports why solvent modeling can be structurally uneven across studies. |
| Surfactants | Indirect | Phillips 2017 | Functional-use QSUR precedent exists; direct QED-style surfactant-likeness score not clearly found in this pass. |
| Human drugs | Baseline only | Bickerton 2012; QED line | QED already covers generic drug-likeness. Our human-drug score should be framed as PubChem-category-specific rather than replacing QED conceptually. |
| Animal drugs | Very limited | Huang 2017 and veterinary-use context | Treat as a novel broad-category extension with only adjacent veterinary drug-likeness precedent. |

## What the literature search means for manuscript framing

1. **Pesticides** already have strong direct score literature.
2. **Fragrances** and **flavor/GRAS chemistry** have strong category-focused chemical-space literature.
3. **Cosmetics, surfactants, solvents, flavorants, and fragrances** also have adjacent functional-use QSUR literature.
4. **Food contact substances** and **animal drugs** look much less developed as continuous category-scoring problems.
5. Therefore the manuscript should not claim that every category has an established prior score. Instead it should say:
   - the study extends the QED paradigm into a mixed landscape,
   - where some categories have direct precedents (especially pesticides),
   - some have strong adjacent chemical-space or functional-use precedents (fragrance, flavor/GRAS, cosmetics, surfactants, solvents),
   - and others remain underdeveloped enough that the present work is gap-filling.

## Minimum reference insertions by manuscript section

### Introduction
- Bickerton 2012
- Tice 2001
- Avram 2014
- Ruddigkeit 2014
- Sprous & Salemme 2007
- Medina-Franco 2012
- Isaacs 2016
- Phillips 2017

### Discussion
- Hao 2011
- Gramatica 1999
- Alder 2016
- Jayaseelan 2012 and/or Brunner 2020

## Sources checked in this pass
- Local PDF: `D:/research/pesticde/research_articles/Pest Management Science - 2001 - Tice - Selecting the right compounds for screening  does Lipinski s Rule of 5 for.pdf`
- PMC / PubMed / DOI reads for Avram 2014, Ruddigkeit 2014, Medina-Franco 2012
- Web search and primary-source retrieval for Isaacs 2016, Phillips 2017, Gramatica 1999, Sprous & Salemme 2007, Alder 2016, Brunner 2020
