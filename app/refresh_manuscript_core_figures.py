from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from app.qed_inspired_validation import CATEGORIES, DISPLAY, MODEL_IDS, ROOT

ANALYSIS = ROOT / "results" / "qed_inspired_analysis"
FIGURES = ROOT / "paper" / "figures"
MODELS = ROOT / "app" / "data" / "models"
REBUILD_SUMMARY = ROOT / "paper" / "final_category_rebuild_summary.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def labels() -> list[str]:
    return [DISPLAY[category] for category in CATEGORIES]


def plot_negative_sets() -> None:
    rows = {row["category"]: row for row in read_csv(REBUILD_SUMMARY)}
    if "flavor_fragrance" not in rows:
        merged_positive = len(read_csv(ROOT / "app" / "output" / "final_category_rebuild" / "inputs" / "flavor_fragrance__positive.csv"))
        merged_negative = len(read_csv(ROOT / "app" / "output" / "final_category_rebuild" / "inputs" / "flavor_fragrance__negative_source.csv"))
        rows["flavor_fragrance"] = {
            "positive_count": str(merged_positive),
            "negative_source_count": str(merged_negative),
            "negative_count_after_tanimoto": str(merged_negative),
        }
    source = [int(rows[category]["negative_source_count"]) for category in CATEGORIES]
    retained = [int(rows[category]["negative_count_after_tanimoto"]) for category in CATEGORIES]
    fractions = [100 * kept / total for kept, total in zip(retained, source)]
    x = np.arange(len(CATEGORIES))
    width = 0.34
    fig, ax = plt.subplots(figsize=(14, 9), constrained_layout=True)
    ax.bar(
        x - width / 2,
        source,
        width,
        label="Cross-category comparison source",
        color="#9ECAE1",
        edgecolor="#333333",
        linewidth=0.6,
    )
    bars = ax.bar(
        x + width / 2,
        retained,
        width,
        label="Retained after Tanimoto < 0.3",
        color="#3182BD",
        edgecolor="#333333",
        linewidth=0.6,
    )
    for bar, fraction in zip(bars, fractions):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(source) * 0.018,
            f"{fraction:.1f}%",
            ha="center",
            va="bottom",
            fontsize=11,
            color="black",
        )
    ax.set_title("Negative-set construction for the retained category scores", fontsize=21, weight="bold", pad=14)
    ax.set_ylabel("Molecule count", fontsize=16)
    ax.set_xlabel("Chemical category", fontsize=16, labelpad=14)
    ax.set_xticks(x, labels(), rotation=45, ha="right", fontsize=13)
    ax.tick_params(axis="y", labelsize=13)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=13, ncol=2)
    fig.savefig(FIGURES / "figure1_final_rebuild_workflow.png", dpi=300, facecolor="white")
    plt.close(fig)


def plot_auc() -> None:
    rows = read_csv(ANALYSIS / "bootstrap_confidence_intervals.csv")
    by_category = {(row["category"], row["metric"]): float(row["estimate"]) for row in rows}
    auc = [by_category[(category, "auc")] for category in CATEGORIES]
    qed = [by_category[(category, "qed_auc")] for category in CATEGORIES]
    x = np.arange(len(CATEGORIES)); width = 0.36
    fig, ax = plt.subplots(figsize=(14, 8), constrained_layout=True)
    ax.bar(x - width / 2, auc, width, label="Category score AUC", color="#3B80C9", edgecolor="#333333", linewidth=0.5)
    ax.bar(x + width / 2, qed, width, label="QED AUC", color="#AFAFAF", edgecolor="#333333", linewidth=0.5)
    ax.set_title("Category-score performance versus QED", fontsize=21, weight="bold", pad=14)
    ax.set_ylabel("Area under the ROC curve", fontsize=16)
    ax.set_xlabel("Chemical category", fontsize=16, labelpad=14)
    ax.set_xticks(x, labels(), rotation=45, ha="right", fontsize=13)
    ax.tick_params(axis="y", labelsize=13)
    ax.set_ylim(0, 1.02); ax.grid(axis="y", alpha=0.25); ax.legend(fontsize=13, ncol=2)
    fig.savefig(FIGURES / "figure2_auc_vs_qed.png", dpi=300, facecolor="white")
    plt.close(fig)


def load_models() -> dict[str, dict]:
    return {
        category: json.loads((MODELS / f"{MODEL_IDS[category]}.json").read_text(encoding="utf-8"))
        for category in CATEGORIES
    }


def pattern_names(config: dict) -> list[str]:
    if config["model_type"] == "han_edc":
        return list(config.get("smarts_patterns", {}))
    return list(config.get("selected_patterns", {}))


def descriptor_names(config: dict) -> list[str]:
    if config["model_type"] == "han_edc":
        return list(config.get("selected_props", config.get("descriptor_ranges", {})))
    return list(config.get("selected_props", []))


def descriptor_weight(config: dict) -> float:
    if config["model_type"] == "han_edc":
        return float(config.get("weights", {}).get("property_score", 0.0))
    return float(config.get("best_w", 0.5))


def plot_composition(models: dict[str, dict]) -> None:
    descriptor_counts = [len(descriptor_names(models[c])) for c in CATEGORIES]
    pattern_counts = [len(pattern_names(models[c])) for c in CATEGORIES]
    weights = [descriptor_weight(models[c]) for c in CATEGORIES]
    x = np.arange(len(CATEGORIES)); width = 0.36
    fig, ax = plt.subplots(figsize=(14, 9), constrained_layout=True)
    ax.bar(x - width / 2, descriptor_counts, width, color="#66BFA5", edgecolor="#333333", label="Descriptor count")
    ax.bar(x + width / 2, pattern_counts, width, color="#FC8D59", edgecolor="#333333", label="Pattern count")
    ax.set_title("Descriptor and structural-pattern composition of category scores", fontsize=21, weight="bold", pad=14)
    ax.set_ylabel("Number of retained terms", fontsize=16)
    ax.set_xlabel("Chemical category", fontsize=16, labelpad=14)
    ax.set_xticks(x, labels(), rotation=45, ha="right", fontsize=13)
    ax.tick_params(axis="y", labelsize=13); ax.grid(axis="y", alpha=0.25)
    for i, value in enumerate(weights):
        ax.text(i, max(descriptor_counts[i], pattern_counts[i]) + 0.3, f"w={value:.2f}", ha="center", va="bottom", fontsize=11)
    ax.legend(fontsize=13, ncol=2)
    fig.savefig(FIGURES / "figure3_descriptor_pattern_composition.png", dpi=300, facecolor="white")
    plt.close(fig)


def plot_patterns(models: dict[str, dict]) -> None:
    all_patterns = sorted(set().union(*(set(pattern_names(models[c])) for c in CATEGORIES)))
    matrix = np.asarray([[name in pattern_names(models[c]) for name in all_patterns] for c in CATEGORIES], dtype=float)
    fig, ax = plt.subplots(figsize=(18, 7.5), constrained_layout=True)
    ax.imshow(matrix, cmap=plt.matplotlib.colors.ListedColormap(["#F4F7F8", "#2C7FB8"]), vmin=0, vmax=1, aspect="auto")
    ax.set_title("Cross-category distribution of retained structural patterns", fontsize=21, weight="bold", pad=14)
    ax.set_xticks(range(len(all_patterns)), [name.replace("_", " ") for name in all_patterns], rotation=45, ha="left", fontsize=10)
    ax.xaxis.tick_top(); ax.set_yticks(range(len(CATEGORIES)), labels(), fontsize=12)
    ax.set_xticks(np.arange(-0.5, len(all_patterns), 1), minor=True); ax.set_yticks(np.arange(-0.5, len(CATEGORIES), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2); ax.tick_params(which="minor", bottom=False, left=False)
    fig.savefig(FIGURES / "figure4_structural_pattern_comparison.png", dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    models = load_models()
    plot_negative_sets(); plot_auc(); plot_composition(models); plot_patterns(models)
    print(FIGURES / "figure1_final_rebuild_workflow.png")
    print(FIGURES / "figure2_auc_vs_qed.png")
    print(FIGURES / "figure3_descriptor_pattern_composition.png")
    print(FIGURES / "figure4_structural_pattern_comparison.png")


if __name__ == "__main__":
    main()
