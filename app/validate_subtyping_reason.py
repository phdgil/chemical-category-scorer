from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

from build_scoring_models import CHOI_PROPERTY_FUNCS, build_choi_model

APP_DIR = Path(__file__).resolve().parent
INPUT_DIR = APP_DIR / "output" / "pubchem_pipeline" / "full_category_decision" / "inputs"
RESULTS_DIR = APP_DIR.parent / "results" / "subtyping_validation"
TMP_DIR = APP_DIR / "output" / "subtyping_validation"
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
CLUSTER_DISTANCE_CUTOFF = 0.45
PROPERTY_KEYS = list(CHOI_PROPERTY_FUNCS.keys())

TARGET_SPECS: dict[str, dict[str, Any]] = {
    "cosmetics": {
        "related": [
            "fragrances",
            "surfactants",
            "food_additives",
            "flavoring_agents",
            "solvents",
            "lipids",
            "plasticizers",
        ],
        "controls": ["pfas", "solvents"],
    },
    "food_contact_substances": {
        "related": [
            "food_additives",
            "solvents",
            "lipids",
            "surfactants",
            "plasticizers",
            "cosmetics",
            "flavoring_agents",
        ],
        "controls": ["pfas", "solvents"],
    },
    "human_drugs": {
        "related": [
            "animal_drugs",
            "endocrine_disruptors",
            "food_additives",
            "flavoring_agents",
            "cosmetics",
            "solvents",
            "pesticides",
        ],
        "controls": ["pfas", "solvents"],
    },
    "animal_drugs": {
        "related": [
            "human_drugs",
            "endocrine_disruptors",
            "food_additives",
            "cosmetics",
            "solvents",
            "pesticides",
            "flavoring_agents",
        ],
        "controls": ["pfas", "solvents"],
    },
}



def read_smiles(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return [row.get("SMILES", "").strip() for row in reader if row.get("SMILES", "").strip()]


def load_category_smiles() -> dict[str, list[str]]:
    categories: dict[str, list[str]] = {}
    for path in sorted(INPUT_DIR.glob("*__positive.csv")):
        slug = path.name.replace("__positive.csv", "")
        categories[slug] = list(dict.fromkeys(read_smiles(path)))
    return categories


def mol_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return mol


def build_entries(smiles_list: list[str]) -> list[dict[str, Any]]:
    entries = []
    for smiles in smiles_list:
        mol = mol_from_smiles(smiles)
        if mol is None:
            continue
        props = {name: float(func(mol)) for name, func in CHOI_PROPERTY_FUNCS.items()}
        entries.append(
            {
                "smiles": smiles,
                "mol": mol,
                "fp": FP_GEN.GetFingerprint(mol),
                "props": props,
            }
        )
    return entries


def write_smiles_csv(path: Path, smiles_list: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["SMILES"])
        for smiles in smiles_list:
            writer.writerow([smiles])


def sample_list(items: list[str], limit: int, seed: int) -> list[str]:
    unique = list(dict.fromkeys(items))
    if len(unique) <= limit:
        return unique
    rng = random.Random(seed)
    return rng.sample(unique, limit)


def property_distance(props: dict[str, float], medians: dict[str, float], scales: dict[str, float]) -> float:
    total = 0.0
    for key in PROPERTY_KEYS:
        scale = max(scales.get(key, 1.0), 1e-6)
        total += abs(props[key] - medians[key]) / scale
    return total / len(PROPERTY_KEYS)


def build_property_matched_negatives(
    positive_entries: list[dict[str, Any]],
    candidate_entries: list[dict[str, Any]],
    limit: int,
) -> list[str]:
    medians = {}
    scales = {}
    for key in PROPERTY_KEYS:
        values = sorted(entry["props"][key] for entry in positive_entries)
        medians[key] = values[len(values) // 2]
        q1 = values[int(0.25 * (len(values) - 1))]
        q3 = values[int(0.75 * (len(values) - 1))]
        scales[key] = max(q3 - q1, 1.0)
    ranked = sorted(
        candidate_entries,
        key=lambda entry: (property_distance(entry["props"], medians, scales), len(entry["smiles"])),
    )
    return [entry["smiles"] for entry in ranked[:limit]]


def butina_clusters(fps: list[Any], cutoff: float = CLUSTER_DISTANCE_CUTOFF) -> list[list[int]]:
    if not fps:
        return []
    if len(fps) == 1:
        return [[0]]
    dists = []
    for idx in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[idx], fps[:idx])
        dists.extend(1.0 - sim for sim in sims)
    return [list(cluster) for cluster in Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)]


def coverage_clusters(cluster_sizes: list[int], fraction: float) -> int:
    total = sum(cluster_sizes)
    running = 0
    for index, size in enumerate(cluster_sizes, start=1):
        running += size
        if running / max(total, 1) >= fraction:
            return index
    return len(cluster_sizes)


def category_structure_metrics(entries: list[dict[str, Any]], seed: int, sample_limit: int = 400) -> dict[str, Any]:
    rng = random.Random(seed)
    sampled = list(entries)
    if len(sampled) > sample_limit:
        sampled = rng.sample(sampled, sample_limit)
    fps = [entry["fp"] for entry in sampled]
    clusters = butina_clusters(fps)
    cluster_sizes = sorted((len(cluster) for cluster in clusters), reverse=True)
    within_nearest = []
    pairwise_means = []
    for idx, fp in enumerate(fps):
        others = fps[:idx] + fps[idx + 1 :]
        if not others:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, others)
        within_nearest.append(max(sims))
        pairwise_means.append(sum(sims) / len(sims))
    return {
        "entry_count": len(entries),
        "sampled_count": len(sampled),
        "cluster_count": len(cluster_sizes),
        "largest_cluster_fraction": round(cluster_sizes[0] / max(len(sampled), 1), 4) if cluster_sizes else 0.0,
        "clusters_for_50pct": coverage_clusters(cluster_sizes, 0.50),
        "clusters_for_80pct": coverage_clusters(cluster_sizes, 0.80),
        "median_nearest_within_similarity": round(sorted(within_nearest)[len(within_nearest) // 2], 4) if within_nearest else 0.0,
        "mean_pairwise_similarity": round(sum(pairwise_means) / len(pairwise_means), 4) if pairwise_means else 0.0,
    }


def best_external_category_distribution(
    positive_entries: list[dict[str, Any]],
    related_entries: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    counts: Counter[str] = Counter()
    best_sims: dict[str, list[float]] = {slug: [] for slug in related_entries}
    for entry in positive_entries:
        best_slug = None
        best_sim = -1.0
        for slug, candidates in related_entries.items():
            if not candidates:
                continue
            sims = DataStructs.BulkTanimotoSimilarity(entry["fp"], [candidate["fp"] for candidate in candidates])
            local_best = max(sims) if sims else 0.0
            if local_best > best_sim:
                best_sim = local_best
                best_slug = slug
        if best_slug is not None:
            counts[best_slug] += 1
            best_sims[best_slug].append(best_sim)
    total = sum(counts.values())
    rows = []
    for slug, count in counts.most_common():
        sims = sorted(best_sims[slug])
        rows.append(
            {
                "best_external_category": slug,
                "count": count,
                "share": round(count / max(total, 1), 4),
                "median_best_similarity": round(sims[len(sims) // 2], 4) if sims else 0.0,
            }
        )
    return rows


def build_validation_model(
    slug: str,
    regime: str,
    positive_smiles: list[str],
    negative_smiles: list[str],
    bayes_trials: int,
    seed: int,
) -> dict[str, Any]:
    run_dir = TMP_DIR / slug / regime
    positive_csv = run_dir / "positive.csv"
    negative_csv = run_dir / "negative.csv"
    model_path = run_dir / "model.json"
    write_smiles_csv(positive_csv, positive_smiles)
    write_smiles_csv(negative_csv, negative_smiles)
    build_choi_model(
        positive_csv=positive_csv,
        negative_source_csv=negative_csv,
        model_id=f"validate_{slug}_{regime}",
        label=f"{slug} {regime}",
        category=slug,
        output_path=model_path,
        bayes_trials=bayes_trials,
        seed=seed,
    )
    return json.loads(model_path.read_text(encoding="utf-8"))


def negative_regimes_for_target(
    slug: str,
    categories: dict[str, list[str]],
    category_entries: dict[str, list[dict[str, Any]]],
    max_negative: int,
    seed: int,
) -> dict[str, list[str]]:
    related = set(TARGET_SPECS[slug]["related"])
    all_other_smiles = []
    for other_slug, smiles in categories.items():
        if other_slug != slug:
            all_other_smiles.extend(smiles)
    all_other_smiles = list(dict.fromkeys(all_other_smiles))

    related_smiles = []
    for other_slug in TARGET_SPECS[slug]["related"]:
        related_smiles.extend(categories.get(other_slug, []))
    related_smiles = list(dict.fromkeys(related_smiles))

    distant_smiles = []
    for other_slug, smiles in categories.items():
        if other_slug == slug or other_slug in related:
            continue
        distant_smiles.extend(smiles)
    distant_smiles = list(dict.fromkeys(distant_smiles))

    candidate_entries = []
    for other_slug, entries in category_entries.items():
        if other_slug != slug:
            candidate_entries.extend(entries)
    property_matched = build_property_matched_negatives(category_entries[slug], candidate_entries, max_negative)

    regimes = {
        "all_other_random": sample_list(all_other_smiles, max_negative, seed + 11),
        "related_hard": sample_list(related_smiles, max_negative, seed + 22),
        "distant_other": sample_list(distant_smiles, max_negative, seed + 33),
        "property_matched": property_matched,
    }
    return regimes


def summarize_target(
    slug: str,
    categories: dict[str, list[str]],
    category_entries: dict[str, list[dict[str, Any]]],
    bayes_trials: int,
    max_negative: int,
    seed: int,
) -> dict[str, Any]:
    structure = category_structure_metrics(category_entries[slug], seed=seed)
    related_entries = {related_slug: category_entries[related_slug] for related_slug in TARGET_SPECS[slug]["related"] if related_slug in category_entries}
    external_distribution = best_external_category_distribution(category_entries[slug], related_entries)
    regimes = negative_regimes_for_target(slug, categories, category_entries, max_negative=max_negative, seed=seed)
    regime_rows = []
    for offset, (regime, negative_smiles) in enumerate(regimes.items(), start=1):
        config = build_validation_model(
            slug=slug,
            regime=regime,
            positive_smiles=categories[slug],
            negative_smiles=negative_smiles,
            bayes_trials=bayes_trials,
            seed=seed + offset,
        )
        metrics = config.get("metrics", {})
        regime_rows.append(
            {
                "category": slug,
                "regime": regime,
                "positive_count": len(categories[slug]),
                "negative_count": int(metrics.get("negative_count", 0)),
                "auc": round(float(metrics.get("auc", 0.0)), 4),
                "ks": round(float(metrics.get("ks", 0.0)), 4),
                "balanced_accuracy": round(float(metrics.get("balanced_accuracy", 0.0)), 4),
                "overlap": round(float(metrics.get("overlap", 0.0)), 4),
                "threshold": round(float(config.get("threshold", 0.5)), 4),
                "selected_props": ",".join(config.get("selected_props", [])),
                "selected_patterns": ",".join(sorted(config.get("selected_patterns", {}).keys())),
                "optimization_method": config.get("optimization_method", ""),
            }
        )
    aucs = [row["auc"] for row in regime_rows]
    summary = {
        "category": slug,
        "structure": structure,
        "external_best_category_distribution": external_distribution,
        "negative_set_validation": regime_rows,
        "auc_range": round(max(aucs) - min(aucs), 4) if aucs else 0.0,
        "best_auc": round(max(aucs), 4) if aucs else 0.0,
        "worst_auc": round(min(aucs), 4) if aucs else 0.0,
    }
    return summary


def compare_controls(
    control_slugs: list[str],
    category_entries: dict[str, list[dict[str, Any]]],
    seed: int,
) -> list[dict[str, Any]]:
    rows = []
    for index, slug in enumerate(control_slugs):
        rows.append({"category": slug, **category_structure_metrics(category_entries[slug], seed=seed + index)})
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate why heterogeneous PubChem categories need subtype evidence rather than one score.")
    parser.add_argument("--bayes-trials", type=int, default=6)
    parser.add_argument("--max-negative", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    categories = load_category_smiles()
    category_entries = {slug: build_entries(smiles) for slug, smiles in categories.items()}

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    target_summaries = []
    negative_rows = []
    analog_rows = []
    for index, slug in enumerate(TARGET_SPECS):
        summary = summarize_target(
            slug=slug,
            categories=categories,
            category_entries=category_entries,
            bayes_trials=args.bayes_trials,
            max_negative=args.max_negative,
            seed=args.seed + 100 * index,
        )
        target_summaries.append(summary)
        negative_rows.extend(summary["negative_set_validation"])
        for row in summary["external_best_category_distribution"]:
            analog_rows.append({"category": slug, **row})

    structure_rows = []
    structure_targets = list(TARGET_SPECS.keys()) + ["pfas", "solvents", "food_additives", "fragrances", "surfactants"]
    for index, slug in enumerate(structure_targets):
        structure_rows.append({"category": slug, **category_structure_metrics(category_entries[slug], seed=args.seed + 500 + index)})

    control_rows = compare_controls(["pfas", "solvents"], category_entries, seed=args.seed + 900)

    summary = {
        "targets": target_summaries,
        "controls": control_rows,
        "interpretation": {
            "cosmetics": "Validate whether cosmetics behaves like a mixture of fragrance/surfactant/solvent/lipid-like families and whether single-score performance is unstable across negative sets.",
            "food_contact_substances": "Validate whether food-contact substances behaves like a broad regulatory bucket with strong overlap to food-additive/solvent/lipid/surfactant families and weak single-score separability.",
            "human_drugs": "Validate whether human drugs behaves like a pharmacology mixture class with overlap to animal drugs and several neighboring functional-use categories.",
            "animal_drugs": "Validate whether animal drugs behaves like a sparse pharmacology/regulatory bucket that should be explained through prototype families rather than one global score.",
        },
    }

    write_csv(
        RESULTS_DIR / "negative_set_sensitivity.csv",
        negative_rows,
        [
            "category",
            "regime",
            "positive_count",
            "negative_count",
            "auc",
            "ks",
            "balanced_accuracy",
            "overlap",
            "threshold",
            "selected_props",
            "selected_patterns",
            "optimization_method",
        ],
    )
    write_csv(
        RESULTS_DIR / "category_structure_metrics.csv",
        structure_rows,
        [
            "category",
            "entry_count",
            "sampled_count",
            "cluster_count",
            "largest_cluster_fraction",
            "clusters_for_50pct",
            "clusters_for_80pct",
            "median_nearest_within_similarity",
            "mean_pairwise_similarity",
        ],
    )
    write_csv(
        RESULTS_DIR / "external_best_category_distribution.csv",
        analog_rows,
        ["category", "best_external_category", "count", "share", "median_best_similarity"],
    )
    (RESULTS_DIR / "subtyping_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(RESULTS_DIR)


if __name__ == "__main__":
    main()
