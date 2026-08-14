from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from rdkit import Chem, DataStructs
from rdkit.Chem import QED, rdFingerprintGenerator
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, roc_auc_score, roc_curve

from app.algorithm_score_engine import score_smiles

ROOT = Path(__file__).resolve().parent.parent
INPUT_DIR = ROOT / "app" / "output" / "final_category_rebuild" / "inputs"
DEFAULT_OUTPUT_DIR = ROOT / "results" / "qed_inspired_analysis"
DEFAULT_FIGURE_DIR = ROOT / "paper" / "figures"
FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

CATEGORIES = (
    "endocrine_disruptors",
    "flavor_fragrance",
    "pesticides",
    "surfactants",
)

DISPLAY = {
    "endocrine_disruptors": "Endocrine\ndisruptors",
    "flavor_fragrance": "Flavor and\nfragrance",
    "pesticides": "Pesticides",
    "surfactants": "Surfactants",
}

MODEL_IDS = {
    "endocrine_disruptors": "han_endocrine_disruptors",
    "flavor_fragrance": "final_flavor_fragrance",
    "pesticides": "final_pesticides",
    "surfactants": "final_surfactants",
}

EXTERNAL = {
    "endocrine_disruptors": (
        "DEDuCT v3 I–III",
        ROOT / "results" / "external_validation" / "deduct_v3_endocrine" / "endocrine_disruptors_true_external_candidates_deduct_I-III.csv",
    ),
    "pesticides": (
        "Health Canada PMRA",
        ROOT / "results" / "external_validation" / "analysis" / "pesticides" / "pesticides_true_external_candidates.csv",
    ),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_smiles(path: Path, columns: Iterable[str] = ("SMILES", "standardized_smiles", "parent_smiles")) -> list[str]:
    rows = read_csv(path)
    values: list[str] = []
    seen: set[str] = set()
    for row in rows:
        value = ""
        for column in columns:
            candidate = str(row.get(column, "")).strip()
            if candidate:
                value = candidate
                break
        if value and value not in seen:
            seen.add(value)
            values.append(value)
    return values


def category_positive_smiles(category: str) -> list[str]:
    return load_smiles(INPUT_DIR / f"{category}__positive.csv")


def score_vector(smiles_values: list[str], model_id: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, list[str]]:
    combined: list[float] = []
    descriptor: list[float] = []
    structural: list[float] = []
    valid_smiles: list[str] = []
    threshold: float | None = None
    for smiles in smiles_values:
        result = score_smiles(smiles, model_id=model_id)
        if not result.valid:
            continue
        valid_smiles.append(smiles)
        combined.append(float(result.score))
        descriptor.append(float(result.property_score))
        structural.append(float(result.structure_score))
        threshold = float(result.threshold)
    if threshold is None:
        raise ValueError(f"No valid scores for {model_id}")
    return (
        np.asarray(combined, dtype=float),
        np.asarray(descriptor, dtype=float),
        np.asarray(structural, dtype=float),
        threshold,
        valid_smiles,
    )


def stable_seed(text: str) -> int:
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16)


def cross_category(output_dir: Path, figure_dir: Path) -> None:
    output_path = output_dir / "cross_category_score_matrix.csv"
    if output_path.is_file():
        rows: list[dict[str, object]] = list(read_csv(output_path))
        rows = [
            row
            for row in rows
            if row["source_category"] in CATEGORIES and row["score_category"] in CATEGORIES
        ]
    else:
        rows = []
        for source_category in CATEGORIES:
            smiles_values = category_positive_smiles(source_category)
            for score_category in CATEGORIES:
                scores, _descriptor, _structural, threshold, valid = score_vector(
                    smiles_values, MODEL_IDS[score_category]
                )
                rows.append(
                    {
                        "source_category": source_category,
                        "score_category": score_category,
                        "model_id": MODEL_IDS[score_category],
                        "input_count": len(smiles_values),
                        "valid_count": len(valid),
                        "median_score": float(np.median(scores)),
                        "mean_score": float(np.mean(scores)),
                        "lower_quartile": float(np.quantile(scores, 0.25)),
                        "upper_quartile": float(np.quantile(scores, 0.75)),
                        "frozen_threshold": threshold,
                        "fraction_at_or_above_threshold": float(np.mean(scores >= threshold)),
                    }
                )
        write_csv(output_path, rows)

    evaluation_path = output_dir / "evaluation_component_scores.csv"
    comparison_fpr: dict[str, float] = {}
    if evaluation_path.is_file():
        counts = {category: [0, 0] for category in CATEGORIES}
        for evaluation_row in read_csv(evaluation_path):
            category = str(evaluation_row["category"])
            if evaluation_row["set"] != "negative" or category not in counts:
                continue
            counts[category][0] += 1
            counts[category][1] += (
                float(evaluation_row["combined_score"]) >= float(evaluation_row["frozen_threshold"])
            )
        comparison_fpr = {
            category: above / total
            for category, (total, above) in counts.items()
            if total
        }
    if len(comparison_fpr) != len(CATEGORIES):
        comparison_fpr = {}
        for category in CATEGORIES:
            negative_scores, _descriptor, _structural, threshold, _valid = score_vector(
                retained_negatives(category, output_dir), MODEL_IDS[category]
            )
            comparison_fpr[category] = float(np.mean(negative_scores >= threshold))

    for row in rows:
        score_category = str(row["score_category"])
        raw_fraction = float(row["fraction_at_or_above_threshold"])
        row["retained_comparison_fraction_at_or_above_threshold"] = comparison_fpr[score_category]
        row["excess_fraction_over_retained_comparison"] = raw_fraction - comparison_fpr[score_category]
    write_csv(output_path, rows)

    lookup = {(str(row["source_category"]), str(row["score_category"])): row for row in rows}
    hit_matrix = np.asarray(
        [
            [float(lookup[(source, score)]["fraction_at_or_above_threshold"]) for score in CATEGORIES]
            for source in CATEGORIES
        ]
    )
    excess_matrix = np.asarray(
        [
            [float(lookup[(source, score)]["excess_fraction_over_retained_comparison"]) for score in CATEGORIES]
            for source in CATEGORIES
        ]
    )

    fig, axes = plt.subplots(1, 2, figsize=(21, 9.5), constrained_layout=True)
    labels = [DISPLAY[category] for category in CATEGORIES]
    limit = max(abs(float(np.min(excess_matrix))), abs(float(np.max(excess_matrix))))
    panels = (
        (
            hit_matrix,
            "Fraction at or above the frozen threshold",
            "magma",
            0.0,
            1.0,
            "Threshold-positive fraction",
            lambda value: f"{100 * value:.0f}%",
        ),
        (
            excess_matrix,
            "Excess over retained-comparison response (percentage points)",
            "RdBu_r",
            -limit,
            limit,
            "Difference from retained comparison (fraction)",
            lambda value: f"{100 * value:+.0f}",
        ),
    )
    for ax, (matrix, title, cmap, vmin, vmax, color_label, formatter) in zip(axes, panels):
        image = ax.imshow(matrix, vmin=vmin, vmax=vmax, cmap=cmap, aspect="auto")
        ax.set_title(title, fontsize=17, weight="bold", pad=16)
        ax.set_xticks(range(len(CATEGORIES)), labels, rotation=45, ha="right", fontsize=11)
        ax.set_yticks(range(len(CATEGORIES)), labels, fontsize=11)
        ax.set_xlabel("Applied category score", fontsize=14, labelpad=10)
        ax.set_ylabel("Source positive set", fontsize=14, labelpad=10)
        for row_index in range(len(CATEGORIES)):
            for column_index in range(len(CATEGORIES)):
                value = matrix[row_index, column_index]
                red, green, blue, _alpha = image.cmap(image.norm(value))
                luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
                color = "black" if luminance > 0.58 else "white"
                ax.text(
                    column_index,
                    row_index,
                    formatter(value),
                    ha="center",
                    va="center",
                    fontsize=9,
                    color=color,
                )
            ax.add_patch(
                Rectangle(
                    (row_index - 0.5, row_index - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor="#00E5FF",
                    linewidth=2.2,
                )
            )
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
        colorbar.set_label(color_label, fontsize=12)
        colorbar.ax.tick_params(labelsize=10)
    fig.suptitle(
        "Threshold response of frozen scores across chemical categories",
        fontsize=21,
        weight="bold",
    )
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / "figure6_cross_category_score_matrix.png", dpi=300, facecolor="white")
    plt.close(fig)
    print(output_path)
    print(figure_dir / "figure6_cross_category_score_matrix.png")


def retained_negatives(category: str, output_dir: Path) -> list[str]:
    cache = output_dir / "retained_negatives" / f"{category}.csv"
    if cache.is_file():
        return load_smiles(cache)
    positives = category_positive_smiles(category)
    negative_source = load_smiles(INPUT_DIR / f"{category}__negative_source.csv")
    positive_fps: list[object] = []
    valid_positives: list[str] = []
    for smiles in positives:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            valid_positives.append(smiles)
            positive_fps.append(FP_GENERATOR.GetFingerprint(mol))
    if not positive_fps:
        raise ValueError(f"No valid positives for {category}")
    threshold_result = score_smiles(valid_positives[0], MODEL_IDS[category])
    model_threshold = float(threshold_result.threshold)
    del model_threshold
    retained: list[str] = []
    for smiles in negative_source:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        fingerprint = FP_GENERATOR.GetFingerprint(mol)
        if max(DataStructs.BulkTanimotoSimilarity(fingerprint, positive_fps)) < 0.3:
            retained.append(smiles)
    write_csv(cache, [{"SMILES": smiles} for smiles in retained])
    return retained


def ensure_evaluation_scores(output_dir: Path) -> Path:
    output_path = output_dir / "evaluation_component_scores.csv"
    if output_path.is_file():
        return output_path
    rows: list[dict[str, object]] = []
    for category in CATEGORIES:
        model_id = MODEL_IDS[category]
        groups = {
            "positive": category_positive_smiles(category),
            "negative": retained_negatives(category, output_dir),
        }
        for group, smiles_values in groups.items():
            combined, descriptor, structural, threshold, valid = score_vector(smiles_values, model_id)
            for smiles, combined_score, descriptor_score, structural_score in zip(
                valid, combined, descriptor, structural
            ):
                mol = Chem.MolFromSmiles(smiles)
                qed_score = float(QED.qed(mol)) if mol is not None else float("nan")
                rows.append(
                    {
                        "category": category,
                        "model_id": model_id,
                        "set": group,
                        "smiles": smiles,
                        "combined_score": float(combined_score),
                        "descriptor_score": float(descriptor_score),
                        "structural_score": float(structural_score),
                        "qed_score": qed_score,
                        "frozen_threshold": threshold,
                    }
                )
    write_csv(output_path, rows)
    return output_path


def component_metrics(labels: np.ndarray, scores: np.ndarray) -> tuple[float, float, float, float]:
    auc = float(roc_auc_score(labels, scores))
    fpr, tpr, thresholds = roc_curve(labels, scores)
    balanced = 0.5 * (tpr + 1.0 - fpr)
    index = int(np.argmax(balanced))
    threshold = float(thresholds[index]) if np.isfinite(thresholds[index]) else 0.5
    predictions = scores >= threshold
    return auc, float(balanced[index]), threshold, float(matthews_corrcoef(labels, predictions))


def ablation(output_dir: Path, figure_dir: Path) -> None:
    score_path = ensure_evaluation_scores(output_dir)
    score_rows = read_csv(score_path)
    rows: list[dict[str, object]] = []
    for category in CATEGORIES:
        category_rows = [row for row in score_rows if row["category"] == category]
        labels = np.asarray([1 if row["set"] == "positive" else 0 for row in category_rows], dtype=int)
        for component, column in (
            ("Descriptor component", "descriptor_score"),
            ("Structural evidence component", "structural_score"),
            ("Combined score", "combined_score"),
        ):
            scores = np.asarray([float(row[column]) for row in category_rows], dtype=float)
            auc, balanced, threshold, mcc = component_metrics(labels, scores)
            rows.append(
                {
                    "category": category,
                    "component": component,
                    "positive_count": int(labels.sum()),
                    "negative_count": int((1 - labels).sum()),
                    "auc": auc,
                    "maximum_balanced_accuracy": balanced,
                    "component_optimal_threshold": threshold,
                    "mcc_at_component_optimal_threshold": mcc,
                }
            )
    output_path = output_dir / "component_ablation.csv"
    write_csv(output_path, rows)

    components = ("Descriptor component", "Structural evidence component", "Combined score")
    colors = ("#4C78A8", "#F58518", "#54A24B")
    width = 0.24
    x = np.arange(len(CATEGORIES))
    fig, ax = plt.subplots(figsize=(18, 8.5), constrained_layout=True)
    for offset, component, color in zip((-width, 0.0, width), components, colors):
        values = [
            float(next(row["auc"] for row in rows if row["category"] == category and row["component"] == component))
            for category in CATEGORIES
        ]
        ax.bar(x + offset, values, width=width, label=component, color=color, edgecolor="#333333", linewidth=0.5)
    ax.axhline(0.5, color="#555555", linestyle="--", linewidth=1.3)
    ax.set_ylim(0.0, 1.03)
    ax.set_xticks(x, [DISPLAY[category] for category in CATEGORIES], rotation=45, ha="right", fontsize=12)
    ax.set_yticks(np.arange(0, 1.01, 0.1))
    ax.tick_params(axis="y", labelsize=12)
    ax.set_xlabel("Chemical category", fontsize=15, labelpad=12)
    ax.set_ylabel("Area under the ROC curve", fontsize=15, labelpad=12)
    ax.set_title("Contribution of descriptor and structural evidence components", fontsize=21, weight="bold", pad=18)
    ax.legend(fontsize=12, ncol=3, loc="upper center")
    ax.grid(axis="y", alpha=0.25)
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / "figure7_component_ablation.png", dpi=300, facecolor="white")
    plt.close(fig)
    print(output_path)
    print(figure_dir / "figure7_component_ablation.png")


def bootstrap(output_dir: Path, replicates: int) -> None:
    score_rows = read_csv(ensure_evaluation_scores(output_dir))
    output: list[dict[str, object]] = []
    for category in CATEGORIES:
        category_rows = [row for row in score_rows if row["category"] == category]
        positive = [row for row in category_rows if row["set"] == "positive"]
        negative = [row for row in category_rows if row["set"] == "negative"]
        positive_scores = np.asarray([float(row["combined_score"]) for row in positive])
        negative_scores = np.asarray([float(row["combined_score"]) for row in negative])
        positive_qed = np.asarray([float(row["qed_score"]) for row in positive])
        negative_qed = np.asarray([float(row["qed_score"]) for row in negative])
        threshold = float(category_rows[0]["frozen_threshold"])
        labels = np.r_[np.ones(len(positive_scores), dtype=int), np.zeros(len(negative_scores), dtype=int)]
        scores = np.r_[positive_scores, negative_scores]
        qed_scores = np.r_[positive_qed, negative_qed]
        predictions = scores >= threshold
        estimates = {
            "auc": float(roc_auc_score(labels, scores)),
            "qed_auc": float(roc_auc_score(labels, qed_scores)),
            "auc_delta_vs_qed": float(roc_auc_score(labels, scores) - roc_auc_score(labels, qed_scores)),
            "balanced_accuracy_at_frozen_threshold": float(balanced_accuracy_score(labels, predictions)),
            "mcc_at_frozen_threshold": float(matthews_corrcoef(labels, predictions)),
        }
        distributions = {key: [] for key in estimates}
        generator = np.random.default_rng(stable_seed(category))
        for _ in range(replicates):
            positive_indices = generator.integers(0, len(positive_scores), len(positive_scores))
            negative_indices = generator.integers(0, len(negative_scores), len(negative_scores))
            boot_scores = np.r_[positive_scores[positive_indices], negative_scores[negative_indices]]
            boot_qed = np.r_[positive_qed[positive_indices], negative_qed[negative_indices]]
            boot_labels = labels
            model_auc = float(roc_auc_score(boot_labels, boot_scores))
            qed_auc = float(roc_auc_score(boot_labels, boot_qed))
            boot_predictions = boot_scores >= threshold
            distributions["auc"].append(model_auc)
            distributions["qed_auc"].append(qed_auc)
            distributions["auc_delta_vs_qed"].append(model_auc - qed_auc)
            distributions["balanced_accuracy_at_frozen_threshold"].append(
                float(balanced_accuracy_score(boot_labels, boot_predictions))
            )
            distributions["mcc_at_frozen_threshold"].append(
                float(matthews_corrcoef(boot_labels, boot_predictions))
            )
        for metric, estimate in estimates.items():
            values = np.asarray(distributions[metric], dtype=float)
            output.append(
                {
                    "category": category,
                    "metric": metric,
                    "estimate": estimate,
                    "ci_lower_95": float(np.quantile(values, 0.025)),
                    "ci_upper_95": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": replicates,
                    "positive_count": len(positive_scores),
                    "negative_count": len(negative_scores),
                    "frozen_threshold": threshold,
                }
            )
    output_path = output_dir / "bootstrap_confidence_intervals.csv"
    write_csv(output_path, output)
    print(output_path)


def external_distributions(output_dir: Path, figure_dir: Path) -> None:
    evaluation_rows = read_csv(ensure_evaluation_scores(output_dir))
    rows: list[dict[str, object]] = []
    plot_data: dict[str, dict[str, np.ndarray]] = {}
    for category, (source, path) in EXTERNAL.items():
        external_smiles = load_smiles(path, ("standardized_smiles", "SMILES", "parent_smiles"))
        external_scores, _descriptor, _structural, threshold, valid = score_vector(
            external_smiles, MODEL_IDS[category]
        )
        internal_positive = np.asarray(
            [
                float(row["combined_score"])
                for row in evaluation_rows
                if row["category"] == category and row["set"] == "positive"
            ]
        )
        retained_negative = np.asarray(
            [
                float(row["combined_score"])
                for row in evaluation_rows
                if row["category"] == category and row["set"] == "negative"
            ]
        )
        plot_data[category] = {
            "PubChem category set": internal_positive,
            "Retained comparison set": retained_negative,
            "Nonoverlapping external set": external_scores,
        }
        for smiles, score in zip(valid, external_scores):
            rows.append(
                {
                    "category": category,
                    "external_source": source,
                    "smiles": smiles,
                    "score": float(score),
                    "frozen_threshold": threshold,
                    "at_or_above_threshold": int(score >= threshold),
                }
            )
    output_path = output_dir / "external_score_distributions.csv"
    write_csv(output_path, rows)

    fig, axes = plt.subplots(
        1,
        len(EXTERNAL),
        figsize=(6 * len(EXTERNAL), 7.5),
        constrained_layout=True,
        sharey=True,
    )
    axes = np.atleast_1d(axes)
    colors = ("#4C78A8", "#B8B8B8", "#F58518")
    for ax, category in zip(axes, EXTERNAL):
        groups = plot_data[category]
        values = list(groups.values())
        violin = ax.violinplot(values, positions=(1, 2, 3), widths=0.78, showmedians=True, showextrema=False)
        for body, color in zip(violin["bodies"], colors):
            body.set_facecolor(color)
            body.set_edgecolor("#333333")
            body.set_alpha(0.78)
        violin["cmedians"].set_color("black")
        violin["cmedians"].set_linewidth(2.0)
        threshold = float(next(row["frozen_threshold"] for row in rows if row["category"] == category))
        ax.axhline(threshold, color="#C00000", linestyle="--", linewidth=1.8, label="Frozen threshold")
        group_labels = (
            f"PubChem\ncategory set\n(n={len(values[0]):,})",
            f"Retained\ncomparison\n(n={len(values[1]):,})",
            f"External\nnonoverlap\n(n={len(values[2]):,})",
        )
        ax.set_xticks((1, 2, 3), group_labels, fontsize=10)
        ax.set_title(DISPLAY[category].replace("\n", " "), fontsize=17, weight="bold", pad=14)
        ax.set_xlabel(f"External source: {EXTERNAL[category][0]}", fontsize=11, labelpad=12)
        ax.grid(axis="y", alpha=0.25)
        ax.set_ylim(-0.02, 1.02)
    axes[0].set_ylabel("Frozen category score", fontsize=15, labelpad=12)
    axes[0].tick_params(axis="y", labelsize=12)
    axes[0].legend(loc="upper right", fontsize=11)
    fig.suptitle("Score distributions for nonoverlapping external positive sets", fontsize=21, weight="bold")
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / "figure8_external_score_distributions.png", dpi=300, facecolor="white")
    plt.close(fig)
    print(output_path)
    print(figure_dir / "figure8_external_score_distributions.png")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        required=True,
        choices=("cross-category", "ablation", "bootstrap", "external"),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.stage == "cross-category":
        cross_category(args.output_dir, args.figure_dir)
    elif args.stage == "ablation":
        ablation(args.output_dir, args.figure_dir)
    elif args.stage == "bootstrap":
        bootstrap(args.output_dir, args.bootstrap_replicates)
    elif args.stage == "external":
        external_distributions(args.output_dir, args.figure_dir)


if __name__ == "__main__":
    main()
