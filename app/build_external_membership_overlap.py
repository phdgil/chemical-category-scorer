from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle


DATASETS = (
    {
        "category": "endocrine_disruptors",
        "display": "Endocrine disruptors",
        "source": "DEDuCT v3 categories I–III",
        "audit": "deduct_v3_endocrine/endocrine_disruptors_external_overlap_audit.csv",
        "identifier_column": "DEDuCT Identifier",
        "pubchem_unique_parents": 5791,
    },
    {
        "category": "pesticides",
        "display": "Pesticides",
        "source": "Health Canada PMRA (resolved structures)",
        "audit": "analysis/pesticides/pesticides_external_overlap_audit.csv",
        "identifier_column": None,
        "pubchem_unique_parents": 2737,
    },
)


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def deduct_i_iii_identifiers(toxicity_path: Path) -> set[str]:
    study_types: dict[str, set[str]] = defaultdict(set)
    for row in read_csv(toxicity_path, delimiter="\t"):
        study_types[row["DEDuCT identifier"]].add(row["Study type"])
    return {
        identifier
        for identifier, types in study_types.items()
        if "IVH" in types or types == {"IVR", "IVTH"} or types == {"IVR"}
    }


def calculate(
    external_root: Path,
    toxicity_path: Path,
) -> list[dict[str, object]]:
    deduct_keep = deduct_i_iii_identifiers(toxicity_path)
    output: list[dict[str, object]] = []
    for dataset in DATASETS:
        rows = read_csv(external_root / str(dataset["audit"]))
        identifier_column = dataset["identifier_column"]
        if identifier_column:
            rows = [row for row in rows if row[identifier_column] in deduct_keep]
        resolved = [row for row in rows if row.get("parent_connectivity_key")]
        external_parents = {row["parent_connectivity_key"] for row in resolved}
        shared_parents = {
            row["parent_connectivity_key"]
            for row in resolved
            if row.get("overlap_scope") == "target_positive"
        }
        pubchem_parents = int(dataset["pubchem_unique_parents"])
        external_unique = len(external_parents)
        shared_unique = len(shared_parents)
        output.append(
            {
                "category": dataset["category"],
                "category_display": dataset["display"],
                "external_source": dataset["source"],
                "external_records": len(rows),
                "structure_resolved_records": len(resolved),
                "external_unique_parent_structures": external_unique,
                "pubchem_unique_parent_structures": pubchem_parents,
                "shared_unique_parent_structures": shared_unique,
                "pubchem_only_unique_parent_structures": pubchem_parents - shared_unique,
                "external_only_unique_parent_structures": external_unique - shared_unique,
                "external_parent_coverage_fraction": shared_unique / external_unique,
            }
        )
    return output


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def draw_venn(ax, row: dict[str, object]) -> None:
    left_color = "#4C78A8"
    right_color = "#F58518"
    left = Circle((0.35, 0.52), 0.27, facecolor=left_color, edgecolor="#1F4E79", alpha=0.50, lw=2.0)
    right = Circle((0.69, 0.52), 0.27, facecolor=right_color, edgecolor="#9C4A00", alpha=0.50, lw=2.0)
    ax.add_patch(left)
    ax.add_patch(right)
    ax.text(
        0.19,
        0.52,
        f"{int(row['pubchem_only_unique_parent_structures']):,}\nPubChem\nonly",
        ha="center",
        va="center",
        fontsize=12,
        weight="bold",
    )
    ax.text(
        0.52,
        0.52,
        f"{int(row['shared_unique_parent_structures']):,}\nshared",
        ha="center",
        va="center",
        fontsize=12,
        weight="bold",
    )
    ax.text(
        0.85,
        0.52,
        f"{int(row['external_only_unique_parent_structures']):,}\nexternal\nonly",
        ha="center",
        va="center",
        fontsize=12,
        weight="bold",
    )
    ax.text(
        0.15,
        0.91,
        f"PubChem category\n(n = {int(row['pubchem_unique_parent_structures']):,})",
        ha="center",
        va="center",
        fontsize=12,
        color="#17365D",
    )
    ax.text(
        0.87,
        0.91,
        f"{row['external_source']}\n(n = {int(row['external_unique_parent_structures']):,})",
        ha="center",
        va="center",
        fontsize=12,
        color="#7F3600",
    )
    coverage = 100 * float(row["external_parent_coverage_fraction"])
    ax.text(
        0.52,
        0.10,
        f"External set represented in PubChem\n{coverage:.1f}%",
        ha="center",
        va="center",
        fontsize=13,
        weight="bold",
    )
    ax.set_title(str(row["category_display"]), fontsize=17, weight="bold", pad=14)
    ax.set_xlim(0, 1.04)
    ax.set_ylim(0, 1.04)
    ax.set_aspect("equal")
    ax.axis("off")


def draw_figure(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(rows), figsize=(18, 7.2))
    if len(rows) == 1:
        axes = [axes]
    for ax, row in zip(axes, rows):
        draw_venn(ax, row)
    fig.subplots_adjust(top=0.78, bottom=0.13, wspace=0.16)
    fig.suptitle(
        "Direct parent-structure overlap between PubChem-derived category sets and external positive sources",
        fontsize=21,
        weight="bold",
        y=0.97,
    )
    fig.text(
        0.5,
        0.015,
        "Counts are unique standardized parent structures. Circle areas are schematic and are not proportional to set size.",
        ha="center",
        fontsize=13,
    )
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, required=True)
    parser.add_argument("--toxicity-file", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-figure", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = calculate(args.external_root, args.toxicity_file)
    write_csv(args.output_csv, rows)
    draw_figure(args.output_figure, rows)
    for row in rows:
        print(
            f"{row['category']}: shared={row['shared_unique_parent_structures']} "
            f"external={row['external_unique_parent_structures']} "
            f"coverage={100 * float(row['external_parent_coverage_fraction']):.1f}%"
        )


if __name__ == "__main__":
    main()
