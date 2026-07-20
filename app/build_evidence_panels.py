from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

from build_scoring_models import CHOI_CANDIDATE_PATTERNS, CHOI_PROPERTY_FUNCS

APP_DIR = Path(__file__).resolve().parent
DATA_DIR = APP_DIR / "data"
EVIDENCE_DIR = DATA_DIR / "evidence_panels"
INPUT_DIR = APP_DIR / "output" / "pubchem_pipeline" / "full_category_decision" / "inputs"

FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
CLUSTER_DISTANCE_CUTOFF = 0.45
MAX_FAMILIES = 6
MAX_MEMBERS_RECORDED = 12
TOP_MOTIFS = 4
MOTIF_MIN_FRACTION = 0.20

PANEL_SPECS = {
    "cosmetics": {
        "positive_csv": INPUT_DIR / "cosmetics__positive.csv",
        "related_models": [
            "choi_fragrance",
            "choi_surfactant",
            "pubchem_solvents",
            "pubchem_flavoring_agents",
            "pubchem_food_additives",
        ],
        "notes": [
            "Cosmetics is a heterogeneous use category, so evidence is organized as prototype families rather than hard subtype predictions.",
            "Top related-model scores are supporting evidence only, not final labels.",
        ],
    },
    "food_contact_substances": {
        "positive_csv": INPUT_DIR / "food_contact_substances__positive.csv",
        "related_models": [
            "pubchem_food_additives",
            "pubchem_solvents",
            "pubchem_lipids",
            "choi_surfactant",
        ],
        "notes": [
            "Food-contact substances is a broad regulatory/use bucket, so evidence is organized as prototype families rather than hard subtype predictions.",
            "Top related-model scores are supporting evidence only, not final labels.",
        ],
    },
    "human_drugs": {
        "positive_csv": INPUT_DIR / "human_drugs__positive.csv",
        "related_models": [
            "pubchem_endocrine_disruptors",
            "pubchem_food_additives",
            "pubchem_flavoring_agents",
            "pubchem_solvents",
            "kim_pesticide",
            "lee_pesticide",
        ],
        "notes": [
            "Human drugs is a pharmacology mixture class with broad scaffold diversity, so evidence is organized as prototype families rather than one hard score.",
            "Top related-model scores are supporting evidence only, not final labels.",
        ],
    },
    "animal_drugs": {
        "positive_csv": INPUT_DIR / "animal_drugs__positive.csv",
        "related_models": [
            "pubchem_endocrine_disruptors",
            "pubchem_food_additives",
            "pubchem_flavoring_agents",
            "pubchem_solvents",
            "kim_pesticide",
            "lee_pesticide",
        ],
        "notes": [
            "Animal drugs is a sparse pharmacology/regulatory bucket, so evidence is organized as prototype families rather than one hard score.",
            "Top related-model scores are supporting evidence only, not final labels.",
        ],
    },
}


def _slug(text: str) -> str:
    chars = []
    for ch in (text or "").strip().lower():
        chars.append(ch if ch.isalnum() else "_")
    slug = "".join(chars)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")


def _read_smiles(path: Path) -> list[str]:
    import csv

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = [row.get("SMILES", "").strip() for row in reader]
    return [row for row in rows if row]


def _prepare_entries(smiles_list: list[str]) -> list[dict[str, Any]]:
    entries = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        entries.append(
            {
                "smiles": smiles,
                "fingerprint": FP_GEN.GetFingerprint(mol),
                "props": {name: float(func(mol)) for name, func in CHOI_PROPERTY_FUNCS.items()},
                "motifs": {
                    name: bool(mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)))
                    for name, smarts in CHOI_CANDIDATE_PATTERNS.items()
                    if Chem.MolFromSmarts(smarts) is not None
                },
            }
        )
    return entries


def _cluster_entries(entries: list[dict[str, Any]]) -> list[list[int]]:
    if len(entries) < 2:
        return [[index] for index in range(len(entries))]
    distances = []
    for i in range(1, len(entries)):
        sims = DataStructs.BulkTanimotoSimilarity(entries[i]["fingerprint"], [entries[j]["fingerprint"] for j in range(i)])
        distances.extend(1.0 - sim for sim in sims)
    clusters = Butina.ClusterData(distances, len(entries), CLUSTER_DISTANCE_CUTOFF, isDistData=True)
    return [list(cluster) for cluster in clusters]


def _cluster_medoid(entries: list[dict[str, Any]], cluster: list[int]) -> int:
    if len(cluster) == 1:
        return cluster[0]
    best_index = cluster[0]
    best_score = -1.0
    fps = [entries[idx]["fingerprint"] for idx in cluster]
    for idx in cluster:
        sims = DataStructs.BulkTanimotoSimilarity(entries[idx]["fingerprint"], fps)
        avg = sum(sims) / len(sims)
        if avg > best_score:
            best_score = avg
            best_index = idx
    return best_index


def _family_name(motifs: list[str], family_index: int) -> str:
    if not motifs:
        return f"prototype_family_{family_index}"
    return "_".join(motifs[:2]) + "_family"


def _build_family(entries: list[dict[str, Any]], cluster: list[int], family_index: int) -> dict[str, Any]:
    members = [entries[idx] for idx in cluster]
    medoid_index = _cluster_medoid(entries, cluster)
    medoid_entry = entries[medoid_index]
    motif_counts = Counter()
    for member in members:
        for name, matched in member["motifs"].items():
            if matched:
                motif_counts[name] += 1
    top_motifs = []
    for name, count in motif_counts.most_common(TOP_MOTIFS):
        fraction = count / len(members)
        if fraction >= MOTIF_MIN_FRACTION:
            top_motifs.append({"name": name, "fraction": round(fraction, 4), "smarts": CHOI_CANDIDATE_PATTERNS[name]})

    prop_medians = {}
    for prop_name in CHOI_PROPERTY_FUNCS:
        values = [member["props"][prop_name] for member in members]
        prop_medians[prop_name] = float(median(values))

    similarities = DataStructs.BulkTanimotoSimilarity(medoid_entry["fingerprint"], [member["fingerprint"] for member in members])
    return {
        "family_id": f"family_{family_index:02d}",
        "family_name": _family_name([item["name"] for item in top_motifs], family_index),
        "member_count": len(members),
        "representative_smiles": medoid_entry["smiles"],
        "representative_similarity_within_family": round(sum(similarities) / len(similarities), 4),
        "top_motifs": top_motifs,
        "property_medians": prop_medians,
        "example_smiles": [member["smiles"] for member in members[:MAX_MEMBERS_RECORDED]],
    }


def build_evidence_panel(panel_id: str, positive_csv: Path, related_models: list[str], notes: list[str]) -> Path:
    smiles_list = _read_smiles(positive_csv)
    entries = _prepare_entries(smiles_list)
    clusters = _cluster_entries(entries)
    clusters.sort(key=len, reverse=True)

    families = []
    for family_index, cluster in enumerate(clusters[:MAX_FAMILIES], start=1):
        families.append(_build_family(entries, cluster, family_index))

    panel = {
        "panel_id": panel_id,
        "positive_csv": str(positive_csv),
        "entry_count": len(entries),
        "related_models": related_models,
        "notes": notes,
        "families": families,
    }
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    output_path = EVIDENCE_DIR / f"{panel_id}.json"
    output_path.write_text(json.dumps(panel, ensure_ascii=False, indent=2), encoding="utf-8")
    return output_path


def main() -> None:
    built = []
    for panel_id, spec in PANEL_SPECS.items():
        built.append(build_evidence_panel(panel_id, spec["positive_csv"], spec["related_models"], spec["notes"]))
    for path in built:
        print(path)


if __name__ == "__main__":
    main()
