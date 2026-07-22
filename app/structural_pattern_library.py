"""Shared structural-pattern libraries used by training and validation."""

from __future__ import annotations


CHOI_CANDIDATE_PATTERNS: dict[str, str] = {
    "halogen": "[Cl,Br,F,I]",
    "chloroaromatic": "[Cl]c",
    "fluoroaromatic": "[F]c",
    "CF3": "C(F)(F)F",
    "triazole": "n1cncn1",
    "pyridine": "c1ccncc1",
    "imidazole": "c1cnc[nH]1",
    "nitro": "[N+](=O)[O-]",
    "carbamate": "OC(=O)N",
    "urea": "NC(=O)N",
    "ester": "C(=O)OC",
    "aromatic_oh": "Oc1ccccc1",
    "aromatic_nh2": "Nc1ccccc1",
    "sulfonate": "S(=O)(=O)",
    "quaternary_n": "[N+]",
    "aldehyde": "[CH1](=O)",
    "long_chain": "CCCCCC",
    "glycol": "OCCO",
    "benzophenone": "c1ccccc1C(=O)c1ccccc1",
    "cinnamate": "OC(=O)/C=C/c1ccccc1",
}
