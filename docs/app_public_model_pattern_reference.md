# Public app model pattern reference

This file lists the structural patterns exposed by each public desktop-app score function and gives a verified probe SMILES for each pattern.

Generated from `app/data/app_release_config.json` and `app/data/models/*.json` on 2026-07-20.

## Public model summary

| Category | App model id | Pattern source | Pattern count |
| --- | --- | --- | ---: |
| Animal Drugs | `final_animal_drugs` | `selected_patterns` | 7 |
| Human Drugs | `final_human_drugs` | `selected_patterns` | 13 |
| Cosmetics | `final_cosmetics` | `selected_patterns` | 7 |
| Endocrine Disruptors | `han_endocrine_disruptors` | `smarts_patterns` | 11 |
| Flavoring Agents | `final_flavoring_agents` | `selected_patterns` | 3 |
| Food Additives | `final_food_additives` | `selected_patterns` | 3 |
| Food Contact Substances | `final_food_contact_substances` | `selected_patterns` | 5 |
| Fragrances | `final_fragrances` | `selected_patterns` | 5 |
| Pesticides | `final_pesticides` | `selected_patterns` | 6 |
| Solvents | `final_solvents` | `selected_patterns` | 1 |
| Surfactants | `final_surfactants` | `selected_patterns` | 4 |

## Pattern details and verified app-test SMILES

### Animal Drugs (`final_animal_drugs`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Imidazole | `c1cnc[nH]1` | Imidazole | `c1cnc[nH]1` | unlikely animal drugs (score 0.070440 vs threshold 0.313836); matched: imidazole |
| Glycol | `OCCO` | Ethylene glycol | `OCCO` | likely animal drugs (score 0.635390 vs threshold 0.313836); matched: glycol |
| Ester | `C(=O)OC` | Methyl acetate | `CC(=O)OC` | unlikely animal drugs (score 0.078445 vs threshold 0.313836); matched: ester |
| Long Chain | `CCCCCC` | Dodecane | `CCCCCCCCCCCC` | unlikely animal drugs (score 0.078445 vs threshold 0.313836); matched: long_chain |
| Halogen | `[Cl,Br,F,I]` | Chloromethane | `CCl` | unlikely animal drugs (score 0.058834 vs threshold 0.313836); matched: halogen |
| Aromatic Oh | `Oc1ccccc1` | Phenol | `Oc1ccccc1` | unlikely animal drugs (score 0.039223 vs threshold 0.313836); matched: aromatic_oh |
| Urea | `NC(=O)N` | Urea | `NC(=O)N` | likely animal drugs (score 0.439223 vs threshold 0.313836); matched: urea |

### Human Drugs (`final_human_drugs`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Fluoroaromatic | `[F]c` | Fluorobenzene | `Fc1ccccc1` | unlikely human drugs (score 0.107267 vs threshold 0.171914); matched: fluoroaromatic; halogen |
| Pyridine | `c1ccncc1` | Pyridine | `c1ccncc1` | unlikely human drugs (score 0.118337 vs threshold 0.171914); matched: pyridine |
| Imidazole | `c1cnc[nH]1` | Imidazole | `c1cnc[nH]1` | unlikely human drugs (score 0.152436 vs threshold 0.171914); matched: imidazole |
| Nitro | `[N+](=O)[O-]` | Nitrobenzene | `O=[N+]([O-])c1ccccc1` | likely human drugs (score 0.185688 vs threshold 0.171914); matched: nitro; aromatic_nh2; quaternary_n |
| Carbamate | `OC(=O)N` | Methyl carbamate | `COC(=O)N` | likely human drugs (score 0.205120 vs threshold 0.171914); matched: carbamate; ester |
| Aromatic Oh | `Oc1ccccc1` | Phenol | `Oc1ccccc1` | likely human drugs (score 0.270940 vs threshold 0.171914); matched: aromatic_oh |
| Aromatic Nh2 | `Nc1ccccc1` | Aniline | `Nc1ccccc1` | likely human drugs (score 0.218278 vs threshold 0.171914); matched: aromatic_nh2 |
| Glycol | `OCCO` | Ethylene glycol | `OCCO` | unlikely human drugs (score 0.143004 vs threshold 0.171914); matched: glycol |
| Long Chain | `CCCCCC` | Dodecane | `CCCCCCCCCCCC` | unlikely human drugs (score 0.087663 vs threshold 0.171914); matched: long_chain |
| Urea | `NC(=O)N` | Urea | `NC(=O)N` | unlikely human drugs (score 0.129546 vs threshold 0.171914); matched: urea |
| Quaternary N | `[N+]` | Tetramethylammonium | `C[N+](C)(C)C` | unlikely human drugs (score 0.015260 vs threshold 0.171914); matched: quaternary_n |
| Ester | `C(=O)OC` | Methyl acetate | `CC(=O)OC` | unlikely human drugs (score 0.015260 vs threshold 0.171914); matched: ester |
| Halogen | `[Cl,Br,F,I]` | Chloromethane | `CCl` | unlikely human drugs (score 0.015260 vs threshold 0.171914); matched: halogen |

### Cosmetics (`final_cosmetics`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Quaternary N | `[N+]` | Tetramethylammonium | `C[N+](C)(C)C` | unlikely cosmetics (score 0.310696 vs threshold 0.358088); matched: quaternary_n |
| Benzophenone | `c1ccccc1C(=O)c1ccccc1` | Benzophenone | `O=C(c1ccccc1)c1ccccc1` | unlikely cosmetics (score 0.178743 vs threshold 0.358088); matched: benzophenone |
| Cinnamate | `OC(=O)/C=C/c1ccccc1` | Methyl cinnamate | `COC(=O)/C=C/c1ccccc1` | likely cosmetics (score 0.438101 vs threshold 0.358088); matched: cinnamate; ester |
| Long Chain | `CCCCCC` | Dodecane | `CCCCCCCCCCCC` | unlikely cosmetics (score 0.348730 vs threshold 0.358088); matched: long_chain |
| Ester | `C(=O)OC` | Methyl acetate | `CC(=O)OC` | unlikely cosmetics (score 0.321858 vs threshold 0.358088); matched: ester |
| Glycol | `OCCO` | Ethylene glycol | `OCCO` | unlikely cosmetics (score 0.268115 vs threshold 0.358088); matched: glycol |
| Sulfonate | `S(=O)(=O)` | Ethanesulfonic acid | `CCS(=O)(=O)O` | unlikely cosmetics (score 0.268115 vs threshold 0.358088); matched: sulfonate |

### Endocrine Disruptors (`han_endocrine_disruptors`)

Pattern source: `smarts_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Phenol | `c[OH]` | Phenol | `Oc1ccccc1` | likely endocrine disruptors (score 0.776920 vs threshold 0.664300); matched: phenol |
| Aniline | `c[NH2,NH1,NH0]` | Aniline | `Nc1ccccc1` | likely endocrine disruptors (score 0.852934 vs threshold 0.664300); matched: aniline |
| Halogenated Aromatic | `c[F,Cl,Br,I]` | Fluorobenzene | `Fc1ccccc1` | likely endocrine disruptors (score 0.671733 vs threshold 0.664300); matched: halogenated_aromatic |
| Nitro Aromatic | `c[N+](=O)[O-]` | Nitrobenzene | `O=[N+]([O-])c1ccccc1` | likely endocrine disruptors (score 0.878551 vs threshold 0.664300); matched: aniline; nitro_aromatic |
| Biphenyl | `c1ccccc1-c1ccccc1` | Biphenyl | `c1ccc(cc1)c1ccccc1` | likely endocrine disruptors (score 0.797370 vs threshold 0.664300); matched: biphenyl |
| Steroid Like | `C1CCC2CCCCC2C1` | Androstane | `C1CC2CCC3C4CCCC4CCC3C2C1` | unlikely endocrine disruptors (score 0.619666 vs threshold 0.664300); matched: steroid_like |
| Amide | `[NX3]C(=O)[CX4]` | Acetamide | `CC(=O)N` | unlikely endocrine disruptors (score 0.000407 vs threshold 0.664300); matched: amide |
| Urea | `[NX3]C(=O)[NX3]` | Urea | `NC(=O)N` | unlikely endocrine disruptors (score 0.000400 vs threshold 0.664300); matched: urea |
| Carbamate | `[NX3]C(=O)O` | Methyl carbamate | `COC(=O)N` | unlikely endocrine disruptors (score 0.000366 vs threshold 0.664300); matched: carbamate; ether |
| Sulfonamide | `S(=O)(=O)N` | Benzenesulfonamide | `NS(=O)(=O)c1ccccc1` | unlikely endocrine disruptors (score 0.647531 vs threshold 0.664300); matched: sulfonamide |
| Ether | `[OD2]([#6])[#6]` | Anisole | `COc1ccccc1` | likely endocrine disruptors (score 0.710429 vs threshold 0.664300); matched: ether |

### Flavoring Agents (`final_flavoring_agents`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Aldehyde | `[CH1](=O)` | Acetaldehyde | `CC=O` | likely flavoring agents (score 0.744038 vs threshold 0.530962); matched: aldehyde |
| Cinnamate | `OC(=O)/C=C/c1ccccc1` | Methyl cinnamate | `COC(=O)/C=C/c1ccccc1` | likely flavoring agents (score 0.622628 vs threshold 0.530962); matched: cinnamate; ester |
| Ester | `C(=O)OC` | Methyl acetate | `CC(=O)OC` | likely flavoring agents (score 0.589769 vs threshold 0.530962); matched: ester |

### Food Additives (`final_food_additives`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Aldehyde | `[CH1](=O)` | Acetaldehyde | `CC=O` | likely food additives (score 0.795113 vs threshold 0.493863); matched: aldehyde |
| Cinnamate | `OC(=O)/C=C/c1ccccc1` | Methyl cinnamate | `COC(=O)/C=C/c1ccccc1` | likely food additives (score 0.807387 vs threshold 0.493863); matched: cinnamate; ester |
| Ester | `C(=O)OC` | Methyl acetate | `CC(=O)OC` | likely food additives (score 0.629631 vs threshold 0.493863); matched: ester |

### Food Contact Substances (`final_food_contact_substances`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Benzophenone | `c1ccccc1C(=O)c1ccccc1` | Benzophenone | `O=C(c1ccccc1)c1ccccc1` | unlikely food contact substances (score 0.010641 vs threshold 0.550000); matched: benzophenone |
| Long Chain | `CCCCCC` | Dodecane | `CCCCCCCCCCCC` | likely food contact substances (score 0.705068 vs threshold 0.550000); matched: long_chain |
| Ester | `C(=O)OC` | Methyl acetate | `CC(=O)OC` | likely food contact substances (score 0.679223 vs threshold 0.550000); matched: ester |
| Glycol | `OCCO` | Ethylene glycol | `OCCO` | likely food contact substances (score 0.627534 vs threshold 0.550000); matched: glycol |
| Sulfonate | `S(=O)(=O)` | Ethanesulfonic acid | `CCS(=O)(=O)O` | likely food contact substances (score 0.627534 vs threshold 0.550000); matched: sulfonate |

### Fragrances (`final_fragrances`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Aldehyde | `[CH1](=O)` | Acetaldehyde | `CC=O` | likely fragrances (score 0.731364 vs threshold 0.517323); matched: aldehyde |
| Cinnamate | `OC(=O)/C=C/c1ccccc1` | Methyl cinnamate | `COC(=O)/C=C/c1ccccc1` | likely fragrances (score 0.667761 vs threshold 0.517323); matched: cinnamate; ester |
| Ester | `C(=O)OC` | Methyl acetate | `CC(=O)OC` | likely fragrances (score 0.650875 vs threshold 0.517323); matched: ester |
| Long Chain | `CCCCCC` | Dodecane | `CCCCCCCCCCCC` | likely fragrances (score 0.517323 vs threshold 0.517323); matched: long_chain |
| Benzophenone | `c1ccccc1C(=O)c1ccccc1` | Benzophenone | `O=C(c1ccccc1)c1ccccc1` | unlikely fragrances (score 0.483552 vs threshold 0.517323); matched: benzophenone |

### Pesticides (`final_pesticides`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Chloroaromatic | `[Cl]c` | Chlorobenzene | `Clc1ccccc1` | unlikely pesticides (score 0.498481 vs threshold 0.498509); matched: chloroaromatic |
| Cf3 | `C(F)(F)F` | Benzotrifluoride | `FC(F)(F)c1ccccc1` | unlikely pesticides (score 0.386365 vs threshold 0.498509); matched: CF3 |
| Triazole | `n1cncn1` | 1,2,4-Triazole | `c1ncn[nH]1` | likely pesticides (score 0.528625 vs threshold 0.498509); matched: triazole |
| Nitro | `[N+](=O)[O-]` | Nitrobenzene | `O=[N+]([O-])c1ccccc1` | likely pesticides (score 0.514150 vs threshold 0.498509); matched: nitro |
| Carbamate | `OC(=O)N` | Methyl carbamate | `COC(=O)N` | likely pesticides (score 0.498509 vs threshold 0.498509); matched: carbamate |
| Urea | `NC(=O)N` | Urea | `NC(=O)N` | likely pesticides (score 0.538155 vs threshold 0.498509); matched: urea |

### Solvents (`final_solvents`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Long Chain | `CCCCCC` | Dodecane | `CCCCCCCCCCCC` | likely solvents (score 0.874479 vs threshold 0.502083); matched: long_chain |

### Surfactants (`final_surfactants`)

Pattern source: `selected_patterns`

| Pattern | SMARTS | Probe molecule | Probe SMILES | Expected app result |
| --- | --- | --- | --- | --- |
| Sulfonate | `S(=O)(=O)` | Ethanesulfonic acid | `CCS(=O)(=O)O` | unlikely surfactants (score 0.506357 vs threshold 0.538281); matched: sulfonate |
| Glycol | `OCCO` | Ethylene glycol | `OCCO` | unlikely surfactants (score 0.459709 vs threshold 0.538281); matched: glycol |
| Long Chain | `CCCCCC` | Dodecane | `CCCCCCCCCCCC` | unlikely surfactants (score 0.338524 vs threshold 0.538281); matched: long_chain |
| Quaternary N | `[N+]` | Tetramethylammonium | `C[N+](C)(C)C` | unlikely surfactants (score 0.245410 vs threshold 0.538281); matched: quaternary_n |
