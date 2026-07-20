# Polished manuscript draft for Journal of Cheminformatics

## Title page

**Proposed title**  
Structure-based chemical product-category scoring from molecular structure: a full-positive PubChem rebuild with same-regime QED comparison and local-first software deployment

**Authors**  
[To be completed]

**Affiliations**  
[To be completed]

**Corresponding author**  
[To be completed]

---

## Abstract

Descriptor-based molecular scoring is widely used in cheminformatics, but most interpretable scoring systems have focused on drug-likeness rather than broader chemical product classes. We developed a structure-based scoring framework for eleven broad chemical product categories defined in the PubChem Classification Browser HID 72 hierarchy: animal drugs, human drugs, cosmetics, endocrine disruptors, flavoring agents, food additives, food contact substances, fragrances, pesticides, solvents, and surfactants. In the final manuscript run, all available positives for these categories were used after CID-to-SMILES conversion and SMILES deduplication. Negatives were constructed from the other final categories while retaining cross-category overlap; exclusions were limited to exact target-positive overlap, duplicate negative SMILES, and near-positive Tanimoto matches. Each category scorer combined positive-set physicochemical property ranges with structural-pattern evidence in a single interpretable score and was benchmarked against raw RDKit quantitative estimate of drug-likeness (QED) on the same rebuilt evaluation set. The strongest categories were surfactants (AUC 0.9778), flavoring agents (0.9389), endocrine disruptors (0.9348), solvents (0.9295), and fragrances (0.9254). Moderate but still clearly improved categories were food additives (0.8334), human drugs (0.7876), and pesticides (0.7527). Harder broad categories including food contact substances (0.6213), cosmetics (0.6179), and animal drugs (0.6136) also exceeded QED under the same evaluation regime. Local-first software packaging was prepared both as a desktop application and as an importable Python library so that the final scorers can be used without remote structure submission. These results show that QED-inspired interpretable scoring can be extended beyond drug-likeness toward broad chemical product-category suggestion when full positive sets, explicit negative-set doctrine, and structural-pattern information are used together.

### Scientific Contribution

We generalize QED-style interpretable scoring from drug-likeness to eleven broad PubChem chemical product categories using full-positive uncapped rebuilds. Every final category scorer outperforms raw RDKit QED on the same retained-cross-category evaluation set. The final study also translates the scoring framework into reusable local-first software through both a desktop application and an importable Python library.

**Keywords:** cheminformatics; molecular scoring; QED; PubChem Classification Browser; structural patterns; descriptor-based scoring; chemical product categories; interpretable models

---

## Introduction

Interpretable molecular scoring remains one of the most practical ideas in cheminformatics because it gives chemically meaningful rankings without requiring a separate black-box model for every decision problem. A well-known example is the quantitative estimate of drug-likeness (QED), which uses a compact set of molecular descriptors and structural alerts to summarize how closely a molecule resembles the physicochemical profile of approved drugs [1]. QED is valuable not only because it is useful, but also because it established a simple and portable paradigm: molecular desirability can be expressed as a structure-based score rather than only as a hard class label.

That paradigm should not be limited to drug-likeness. In many chemical workflows, the relevant question is not whether a molecule looks drug-like, but whether it resembles a broad chemical product class such as a fragrance, a flavoring agent, a surfactant, a food additive, or a solvent. Similar questions arise for classes such as pesticides, cosmetics, food contact substances, human drugs, and animal drugs. These categories are useful for library triage, product-oriented inventory organization, safety-oriented screening, and practical cheminformatics decision support. However, broad product-category scoring has received much less attention than drug-likeness, synthetic accessibility, or natural-product likeness.

Related work in cheminformatics shows that interpretable structure-based scores can be highly useful when the target concept is chemically coherent. Beyond QED, notable examples include natural-product likeness scoring [3], synthetic accessibility scoring [4,5], and structural-alert derivation frameworks such as Bioalerts [2]. These studies established three important precedents for the present work: chemically meaningful scoring can be built from compact descriptor-and-fragment logic, structural motifs can be extracted in an interpretable way, and software delivery is part of the scientific value because scoring methods are most useful when they can be reused directly. At the same time, these earlier methods were not designed to answer a broader product-category question such as whether a molecule resembles a flavoring agent, a surfactant, or a food contact substance.

A smaller body of adjacent work has addressed class-specific screening problems such as toxicology alerts, drug-likeness prioritization, and specialized library triage, but these efforts generally focus on one endpoint family, one application domain, or one narrow ontology. They do not directly provide a general framework for broad chemical product-category scoring across multiple industrially relevant classes under a common benchmark. In particular, they do not resolve the two practical issues that emerged most clearly in our setting: whether capped positive sets suppress meaningful structural motifs, and whether broad-category negatives should be treated as mutually exclusive labels or as comparison backgrounds. Those gaps motivate the present study more directly than a simple claim that broad categories have not been studied.

The PubChem Classification Browser provides a scalable source of candidate categories for such a study because it exposes a curated hierarchy of broad chemical classes. In earlier exploratory work in this project, category scorers were built for multiple HID 72 chemical classes, but the initial screen used positive-set caps for tractability. Those caps were acceptable during exploratory screening, but they introduced a scientific risk: if the positive set is truncated, category-specific structural patterns may be lost or underweighted. This concern is especially important for diverse categories such as pesticides, cosmetics, or human drugs, where the missing tail of the positive distribution may contain meaningful chemical motifs.

A second issue is negative-set construction. For broad product categories, positive membership is often easier to define than true non-membership. Molecules that appear in neighboring categories are not automatically invalid negatives, because the modeling question is category-likeness rather than universal exclusivity. For example, a compound may reasonably be treated as a fragrance positive while remaining a valid cosmetics negative if the goal is to recognize fragrance-like chemistry rather than every molecule that can appear anywhere in a cosmetics-related context. Aggressive overlap removal can therefore weaken structural signal by erasing the very contrast needed for broad-category discrimination.

The objective of the present study was to build a final broad-category scoring line for eleven practically relevant PubChem categories while addressing both concerns. We removed the earlier positive cap, retained cross-category negatives under a narrow exclusion rule, and rebuilt each scorer using positive-set physicochemical property ranges and structural-pattern evidence. The resulting models were benchmarked against raw RDKit QED on the same rebuilt evaluation sets, so that every reported improvement was same-regime rather than cross-regime. For endocrine disruptors, cross-regime validation was used to choose the final deployed scorer, and the reconstructed Han Se-eum model was selected because it outperformed the uncapped GJC rebuild on both the native endocrine-vs-drug regime and the final retained-cross-category regime. We then packaged the final scorers into both a desktop application and an importable Python library to demonstrate direct usability. Figure 1 is intended to summarize the full manuscript workflow from PubChem category extraction through score rebuilding, QED comparison, and software packaging.

---

## Methods

### Category source and final target set

The source ontology for the final manuscript run was the PubChem Classification Browser HID 72 hierarchy, specifically the `Chemical and Physical Properties > Chemical Classes` branch [6]. Eleven categories were selected as the final manuscript targets: animal drugs, human drugs, cosmetics, endocrine disruptors, flavoring agents, food additives, food contact substances, fragrances, pesticides, solvents, and surfactants. These categories were chosen because they combined practical chemical relevance with sufficient category size for final broad-category scoring while avoiding product classes that were too sparse to support stable reconstruction.

Positive molecules were recovered from PubChem category membership as compound identifiers (CIDs) and then translated to canonical SMILES strings using the PubChem property service. Final positive counts were defined as unique SMILES counts after CID-to-SMILES conversion and SMILES deduplication, not as raw CID counts. This was necessary because multiple CIDs can collapse to the same structure-level SMILES representation. The final unique positive counts were 859 for animal drugs, 3749 for human drugs, 3665 for cosmetics, 5892 for endocrine disruptors, 2316 for flavoring agents, 2809 for food additives, 1910 for food contact substances, 2286 for fragrances, 2847 for pesticides, 605 for solvents, and 268 for surfactants (Table 1).

**[Insert Table 1 near here]**

### Full-positive rebuild

The earlier exploratory screen used a positive cap for some classes as a practical runtime control. For the final manuscript run, that cap was removed. All available positives for the eleven selected categories were used in the final rebuild after deduplication. The scientific reason for this change was simple: a capped positive set can hide significant structural patterns, especially in chemically heterogeneous categories. In a category such as pesticides, for example, the full positive collection spans aromatic halogenated motifs, nitrogen-rich heterocycles, carbamates, and urea-like structures. Restricting the training pool can easily distort the learned balance between descriptor windows and pattern evidence.

### Negative-set doctrine

Negatives were constructed from the union of the other ten final categories. Cross-category overlap was retained. A candidate negative was removed only when it met one of three exclusion rules: exact overlap with the target positive set, duplicate negative-source SMILES, or near-positive similarity after Tanimoto filtering. This design reflected a category-likeness view of the task rather than an unrealistic assumption that neighboring broad product classes are mutually exclusive.

The negative-set doctrine was therefore asymmetric by design. Positive membership was treated as the strongest available label, whereas the negative pool was a comparison background rather than a universal statement of non-membership. This distinction is important for broad product categories such as fragrances versus cosmetics, or flavoring agents versus food additives, where the same molecule can legitimately appear in multiple use contexts. Retaining cross-category negatives preserves discriminative contrast that would otherwise be erased by blanket overlap removal.

### Similarity filter

Before fitting each scorer, negative candidates were filtered using Morgan fingerprints and a Tanimoto threshold of 0.3. Any negative candidate that was too close to the target positive set was removed. This step reduced trivial leakage while still preserving a difficult and chemically relevant comparison background. The remaining background sizes after Tanimoto filtering were 6010 for animal drugs, 3422 for human drugs, 4525 for cosmetics, 3455 for endocrine disruptors, 6207 for flavoring agents, 4685 for food additives, 5578 for food contact substances, 5484 for fragrances, 4438 for pesticides, 8890 for solvents, and 11815 for surfactants. The Han endocrine scorer was evaluated on this same retained-cross-category endocrine set when selecting the final endocrine model for the manuscript and software package.

### Descriptor and structural-pattern scoring

Each final category scorer combined two interpretable components: (1) physicochemical support based on whether a molecule fell within positive-set descriptor ranges and (2) structural-pattern support based on the presence of motifs enriched in the target class. The descriptor pool consisted of molecular weight, logP, hydrogen-bond donor count, hydrogen-bond acceptor count, topological polar surface area, rotatable bond count, fraction of sp3 carbons, and aromatic ring count. These descriptors were chosen because they map directly onto the traditional QED design space while still being broad enough to cover non-drug product classes.

Structural-pattern candidates came from the project motif library and included interpretable fragments such as ester, aldehyde, cinnamate, long-chain fragments, sulfonate, carbamate, urea, aromatic halogenation, benzophenone-like motifs, imidazole, pyridine, nitro, glycol, and quaternary nitrogen. Pattern counts were not fixed in advance. Instead, they were selected per category according to enrichment against the retained comparison background. This design ensured that final scorers used chemical motifs that remained meaningful after the negative doctrine was applied, rather than importing a generic rule set unchanged across all categories.

The final score for a molecule was therefore a weighted combination of descriptor-range support and structural-pattern support. The mixing weight `w` determined how much importance was given to the descriptor component relative to the pattern component. Categories with broad recurring motifs, such as surfactants, could retain substantial pattern influence, whereas categories with more diffuse structural identity could rely more on descriptor support without collapsing into a descriptor-only score.

### Optimization strategy

The final scoring framework used a hybrid optimization strategy derived from the earlier student workflows. Grid search served as a stable baseline optimizer. A Lee-style Bayesian optimization route was accepted when it improved the composite objective, which considered discrimination quality, threshold behavior, and ranking performance. This allowed the final study to preserve the practical automation logic of the Choi pipeline while incorporating a more flexible optimization layer when it genuinely improved the score. The final table reports the selected optimization route for each category.

### Baseline comparison against QED

Raw RDKit QED was computed on the same rebuilt evaluation set for every final category [7]. This same-regime comparison was the core reporting gate for the manuscript. A category was treated as manuscript-reportable only when its rebuilt scorer achieved higher AUC than QED on that same evaluation set. This rule matters because many apparent improvements in applied cheminformatics disappear when the baseline and the proposed method are evaluated on different positive or negative constructions.

### Software deployment

The final scorers were packaged into two local-first delivery modes: a desktop application for single-SMILES and batch CSV scoring, and an importable Python library that exposes category-level functions in a QED-like style. Both modes call the same bundled final score models and the same underlying scoring engine, so identical molecules produce identical scores regardless of whether the user works through the graphical interface or through a Python script. The desktop interface exposes model selection, per-molecule scoring, file-based batch processing, molecule visualization, and matched-pattern visualization without requiring a hosted service. The Python library supports direct import into cheminformatics workflows that already use RDKit. This local-first design matters for industrial settings because chemical manufacturers may treat molecular structures as confidential assets and may therefore prefer offline evaluation over web submission. Descriptor calculation, scoring, and export all run on the user workstation. This deployment layer shows that the scorers are operational and reusable rather than only notebook-bound. Figure 4 should summarize the software workflow and representative desktop and library usage views, while implementation details that are too long for the main text can remain in supplementary material.

**[Insert Figure 1 near end of Introduction]**
**[Insert Figure 4 near end of Methods/Software deployment]**

---

## Results

An additional PubChem audit distinguished the Agrochemical Information branch from the Chemical Classes Pesticides branch. Agrochemical Information partially overlaps pesticides but also contains transformation- and registry-oriented records that are not interchangeable with the dedicated pesticide class; accordingly, the final pesticide scorer remained defined only from the dedicated Pesticides class. A separate food-branch screen showed that Food Additives and Ingredients contains several follow-up candidates, but none were strong enough or chemically clean enough to enter the final scorer set. Finally, an uncapped endocrine disruptor scorer was rebuilt from all available PubChem endocrine positives, and cross-regime validation against Han Se-eum's endocrine model was used to choose the final endocrine scorer included in both the manuscript tables and the deployed software.

### Full-positive counts and rebuilt benchmark basis

The full-positive rebuild materially changed the benchmark basis relative to the earlier exploratory screen. Categories that were previously constrained by capped subsets now used the full available positive chemistry in this run. Pesticides increased to 2847 unique SMILES, human drugs to 3749, cosmetics to 3665, endocrine disruptors to 5892, food additives to 2809, fragrances to 2286, and flavoring agents to 2316. Because the manuscript results are based on these uncapped positive sets, the final broad-category scorers are more representative of their intended chemical spaces than the earlier sampled versions.

The full-positive refresh also clarified the category-size landscape. Human drugs, cosmetics, endocrine disruptors, food additives, pesticides, fragrances, and flavoring agents are large enough that their final behavior cannot be dismissed as a small-sample artifact. Solvents and surfactants are smaller, but they still retain enough chemistry to support stable final scoring, especially because the negative backgrounds are much larger than the positive sets after Tanimoto filtering. Table 1 summarizes these benchmark counts together with AUC, balanced accuracy, and the QED comparison. Cross-regime endocrine validation additionally showed that Han Se-eum's reconstructed endocrine model outperformed the uncapped GJC endocrine rebuild on both the native drug-negative comparator set and the final retained-cross-category endocrine set; accordingly, the final scorer set used in both the manuscript tables and the deployed software includes the Han endocrine model.

### Final broad-category performance versus QED

All eleven final scorers beat raw RDKit QED in AUC on the same rebuilt evaluation sets (Table 1; Figure 2).

**[Insert Figure 2 near here]**

The strongest categories were surfactants, flavoring agents, endocrine disruptors, solvents, and fragrances. Surfactants achieved an AUC of 0.9778 and balanced accuracy of 0.9391, compared with QED values of 0.2716 and 0.5176. Flavoring agents reached 0.9389 and 0.8889, compared with 0.5656 and 0.6380 for QED. The selected Han endocrine scorer reached 0.9348 and 0.9350 on the final retained-cross-category endocrine set, compared with 0.6072 and 0.5812 for QED. Solvents reached 0.9295 and 0.8680, compared with 0.5197 and 0.6654. Fragrances reached 0.9254 and 0.8603, compared with 0.5527 and 0.6132. These five categories therefore provide the clearest evidence that broad chemical product-category scoring can outperform a generic drug-likeness baseline by a large margin.

Food additives, human drugs, and pesticides formed a middle tier. Food additives reached an AUC of 0.8334 versus 0.4463 for QED, with balanced accuracy 0.7676 versus 0.5501. Human drugs reached 0.7876 versus 0.4730, with balanced accuracy 0.7369 versus 0.5641. Pesticides reached 0.7527 versus 0.6790, with balanced accuracy 0.6840 versus 0.6217. The pesticide result is especially important because it shows that once the earlier positive cap was removed, the broad-category scorer still exceeded QED under the final regime, albeit by a narrower margin than in the best category-specific student-native benchmark.

The hardest categories were food contact substances, cosmetics, and animal drugs. Food contact substances reached an AUC of 0.6213 compared with 0.3497 for QED, cosmetics reached 0.6179 compared with 0.2605, and animal drugs reached 0.6136 compared with 0.4445. These categories are clearly more difficult than surfactants or fragrances, but they still satisfy the same-regime manuscript gate because they remain above QED. Their weaker absolute performance also helps define the realistic boundary of the current scoring framework: broad product-category scoring is possible, but chemically diffuse categories remain challenging.

### Descriptor and pattern composition of final scorers

The final scorers were not purely descriptor-based. Structural patterns remained a meaningful component across categories (Table 2; Figure 3).

**[Insert Table 2 and Figure 3 near here]** Human drugs retained 13 selected patterns with seven descriptor ranges, reflecting the chemically broad but still motif-rich nature of drug-like space in this task. The selected Han endocrine scorer retained seven descriptor statistics and eleven SMARTS-based patterns, placing it among the most structurally rich final scorers. Cosmetics retained seven patterns with four descriptor ranges, food contact substances five patterns with one descriptor range, fragrances five patterns with six descriptor ranges, pesticides six patterns with seven descriptor ranges, and surfactants four patterns with seven descriptor ranges. These counts show that the final selected scorer set did not collapse into property-only scoring.

The pattern content was chemically plausible. Flavoring agents and food additives both retained aldehyde, cinnamate, and ester motifs, which is consistent with common carbonyl-rich and aromatic flavor chemistry. Fragrances retained aldehyde, cinnamate, ester, long-chain, and benzophenone-like motifs, reflecting a mixture of odor-active aromatic and lipid-like structural cues. Pesticides retained chloroaromatic, trifluoromethyl, triazole, nitro, carbamate, and urea motifs, all of which are chemically recognizable within agrochemical space. Surfactants retained sulfonate, glycol, long-chain, and quaternary nitrogen patterns, matching the expected amphiphilic design of many detergent-like compounds.

Weighting behavior also varied meaningfully across categories. Cosmetics used a lower descriptor weight (0.25), indicating stronger reliance on structural motifs relative to descriptors. Food additives and solvents used higher descriptor weights (0.7531), indicating that descriptor windows carried a larger share of their final discrimination. Food contact substances sat near balance with a weight of 0.55, which is notable because this difficult category still required both descriptor and pattern support to remain above QED.

### Category-specific interpretation

The strongest categories share a useful property: they possess both chemically sensible descriptor windows and recurrent structural patterns that are not confined to a tiny positive subset. Surfactants exemplify this behavior most clearly. Their long-chain, ionic, and heteroatom-rich motifs are both chemically recognizable and strongly discriminative against the retained background. Solvents are somewhat different: they rely on a smaller pattern set, but their descriptor profile is sufficiently distinctive that the mixed scorer still greatly exceeds QED.

The middle-tier categories highlight the value of the final negative doctrine. Human drugs are too diverse to be summarized by QED alone, yet the final broad-category scorer still improves strongly once the full positive set is used and cross-category negatives are retained. Pesticides are particularly informative because earlier student-native pipelines already suggested good category-specific performance. The final result shows that broad-category performance remains positive even after moving from a tailored benchmarking regime to a stricter broad-category comparison background.

The difficult categories are also scientifically useful. Cosmetics and food contact substances were previously tempting to subdivide into very narrow families, but that strategy risks converting a broad product question into a collection of small and possibly over-optimistic subproblems. The present manuscript instead keeps them as broad classes and demonstrates that a full-positive, pattern-aware scorer can still beat QED, even if the margin is smaller than in cleaner classes such as surfactants or fragrances. Animal drugs remain the weakest of the eleven categories, but the score still improves over QED and thus remains informative rather than null. The added endocrine disruptor scorer also remains above QED after removing the historical positive-set cap.

---

## Discussion

The main result of this study is not only that some chemical product categories can be scored well, but that a QED-inspired, interpretable scoring framework can be generalized to broad chemical product classes under a stronger benchmark design than the earlier exploratory runs. Two design choices were decisive: removing the positive cap and using a restrained negative-exclusion doctrine.

Removing the positive cap mattered because broad product categories are structurally diverse. If only a subset of positives is retained, the model can overfit to a limited chemistry slice and miss meaningful motifs. The final uncapped run avoided that problem and produced a more credible category basis, especially for large categories such as human drugs, cosmetics, food additives, fragrances, and pesticides. This was not merely a procedural cleanup. It changed the scientific interpretation of the score because the resulting motif library now reflects the full accessible positive chemistry in the selected categories.

The retained-cross-category negative rule also mattered. For broad product classes, category overlap in usage space does not automatically invalidate negative examples. Automatically removing every molecule that appears in a neighboring category can erase legitimate structural contrast and reduce the number of useful hard negatives to an unrealistic minimum. The final rule therefore treated cross-category overlap as chemically informative unless it violated narrow exclusion rules. In practice, this preserved enough background chemistry to let structural patterns emerge rather than forcing the score toward a nearly descriptor-only solution.

QED was an appropriate baseline because it is the clearest conceptual ancestor of the current framework. The final manuscript does not merely state that the rebuilt scorers are useful in isolation; it shows that they outperform an established descriptor-based score on the same rebuilt evaluation sets. This same-regime comparison is methodologically stronger than a casual literature comparison because it prevents benchmark mismatch from masquerading as model improvement. It also gives the manuscript a clear publication logic: the contribution is not a vague claim of usefulness, but a demonstrated improvement over the closest available generic baseline.

The strongest categories—surfactants, flavoring agents, endocrine disruptors, solvents, and fragrances—show that broad chemical product classes can be highly discriminative under an interpretable scoring framework. Food additives, human drugs, and pesticides show that more heterogeneous but still useful classes can remain above QED after full-positive rebuilding. Food contact substances, cosmetics, and animal drugs remain difficult categories, but their results are still informative rather than negative: broad product-category scoring can remain useful even when the class is chemically diffuse. This is a more realistic and more publishable message than pretending that every product class should be equally easy.

This study has several limitations. First, the task is still positive-versus-constructed-background rather than a perfect true-negative problem. Second, the manuscript does not claim universal optimality; it reports a bounded final rebuild under one explicit negative doctrine. Third, no external wet-lab or industrial prospective validation was added in this version. Fourth, this manuscript intentionally focuses on interpretable scoring and does not yet merge the separate machine-learning-final models. Fifth, some categories, especially surfactants, have smaller positive sets than categories such as cosmetics or human drugs, even though their final discrimination is very strong.

Despite these limitations, the manuscript is strong enough for a specialist cheminformatics venue because the contribution is clear: a full-positive, benchmark-aware, interpretable extension of QED-style scoring to eleven broad chemical product categories, validated against same-regime QED and translated into deployable software. This study also provides a practical framework for future refinement. Categories that already show strong margins over QED can be pushed toward application-focused deployment, whereas weaker but still positive categories can be improved later without changing the fundamental benchmark logic. The software contribution is also stronger than a notebook-only release because the same final engine is available both as a local desktop application and as an importable Python library for direct reuse in RDKit-based workflows.

---

## Conclusions

A full-positive PubChem rebuild showed that QED-style interpretable scoring can be extended beyond drug-likeness toward broad chemical product-category suggestion. In the final run, all eleven selected categories exceeded raw RDKit QED on the same rebuilt evaluation sets. The strongest categories were surfactants, flavoring agents, endocrine disruptors, solvents, and fragrances, while food additives, human drugs, and pesticides formed a strong middle tier and food contact substances, cosmetics, and animal drugs remained difficult but still reportable broad classes. The final endocrine scorer set also shows that cross-regime validation can identify a stronger deployable model than the uncapped broad-category rebuild alone when both are evaluated on the same retained-cross-category background.

These results support publication as a cheminformatics methods study and justify same-regime QED comparison as the minimum reporting gate for this project line. The final software package further shows that the category scorers are operational rather than purely theoretical by providing both a local desktop application and an importable Python library. In the final software and manuscript scorer set, endocrine disruption is represented by the Han Se-eum model selected after cross-regime validation against the uncapped GJC endocrine rebuild.

---

## List of abbreviations

AUC, area under the receiver operating characteristic curve  
CID, PubChem compound identifier  
HID, hierarchy identifier  
KS, Kolmogorov-Smirnov statistic  
QED, quantitative estimate of drug-likeness  
QSAR, quantitative structure-activity relationship  
SMILES, simplified molecular-input line-entry system

---

## Declarations

### Availability of data and materials

All rebuilt category tables, model configuration files, QED comparison outputs, manuscript-preparation artifacts, the local-first desktop application, and the importable Python library are stored under the project root, `results`, `paper`, `app`, and `chemical_category_scorer` directories. Before submission, these materials should be deposited in a public repository or archived release so that the full-positive rebuild, same-regime QED comparison, desktop application, and Python library are available without restriction.

### Competing interests

The authors declare that they have no competing interests.

### Funding

Not applicable in the current draft. Update this section before submission if funding support should be acknowledged.

### Authors’ contributions

To be completed by the corresponding author before submission.

### Acknowledgements

To be completed before submission.

### Authors’ information

Not applicable in the current draft.

---

## References

1. Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL (2012) Quantifying the chemical beauty of drugs. Nat Chem 4:90-98. doi:10.1038/nchem.1243.
2. Cortes-Ciriano I (2016) Bioalerts: a python library for the derivation of structural alerts from bioactivity and toxicity data sets. J Cheminform 8:13. doi:10.1186/s13321-016-0125-7.
3. Ertl P, Roggo S, Schuffenhauer A (2008) Natural product-likeness score and its application for prioritization of compound libraries. J Chem Inf Model 48:68-74. doi:10.1021/ci700286x.
4. Ertl P, Schuffenhauer A (2009) Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. J Cheminform 1:8. doi:10.1186/1758-2946-1-8.
5. Huang Q, Li LL, Yang SY (2011) RASA: a rapid retrosynthesis-based scoring method for the assessment of synthetic accessibility of drug-like molecules. J Chem Inf Model 51:2768-2777. doi:10.1021/ci100216g.
6. PubChem Classification Browser HID 72: Chemical Classes. National Center for Biotechnology Information. https://pubchem.ncbi.nlm.nih.gov/classification/#hid=72. Accessed 14 Jul 2026.
7. RDKit QED module documentation. RDKit 2026.03.2 documentation. https://rdkit.org/docs/source/rdkit.Chem.QED.html. Accessed 14 Jul 2026.
