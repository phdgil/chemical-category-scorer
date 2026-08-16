from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt
from rdkit import Chem

from app.qed_inspired_validation import ROOT, category_positive_smiles

ANALYSIS = ROOT / "results" / "qed_inspired_analysis"
PAPER = ROOT / "paper"

DISPLAY = {
    "animal_drugs": "Animal drugs",
    "human_drugs": "Human drugs",
    "cosmetics": "Cosmetics",
    "food_contact_substances": "Food contact",
    "food_additives": "Food additives",
    "solvents": "Solvents",
    "flavoring_agents": "Flavoring agents",
    "fragrances": "Fragrances",
    "endocrine_disruptors": "Endocrine disruptors",
    "pesticides": "Pesticides",
    "surfactants": "Surfactants",
}
CATEGORIES = tuple(DISPLAY)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def canonical_positive_sets() -> dict[str, set[str]]:
    output: dict[str, set[str]] = {}
    for category in CATEGORIES:
        values: set[str] = set()
        for smiles in category_positive_smiles(category):
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is not None:
                values.add(Chem.MolToSmiles(molecule, canonical=True))
        output[category] = values
    return output


def build_overlap_assets() -> tuple[Path, Path, Path, Path]:
    sets = canonical_positive_sets()
    pairwise_rows: list[dict[str, object]] = []
    directional = np.zeros((len(CATEGORIES), len(CATEGORIES)), dtype=float)
    for row_index, category_a in enumerate(CATEGORIES):
        for column_index, category_b in enumerate(CATEGORIES):
            shared = len(sets[category_a] & sets[category_b])
            directional[row_index, column_index] = shared / len(sets[category_a])
            if row_index >= column_index:
                continue
            union = len(sets[category_a] | sets[category_b])
            pairwise_rows.append(
                {
                    "category_a": category_a,
                    "category_b": category_b,
                    "category_a_count": len(sets[category_a]),
                    "category_b_count": len(sets[category_b]),
                    "shared_structure_count": shared,
                    "fraction_of_category_a": shared / len(sets[category_a]),
                    "fraction_of_category_b": shared / len(sets[category_b]),
                    "jaccard_similarity": shared / union,
                }
            )
    pairwise_path = PAPER / "supporting_information_pairwise_category_overlap.csv"
    with pairwise_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pairwise_rows[0]))
        writer.writeheader()
        writer.writerows(pairwise_rows)

    labels = [DISPLAY[category] for category in CATEGORIES]
    fig, ax = plt.subplots(figsize=(14, 11), constrained_layout=True)
    image = ax.imshow(directional, cmap="YlOrRd", vmin=0, vmax=1)
    for row_index in range(len(CATEGORIES)):
        for column_index in range(len(CATEGORIES)):
            value = directional[row_index, column_index]
            ax.text(
                column_index,
                row_index,
                f"{100 * value:.0f}%",
                ha="center",
                va="center",
                fontsize=9,
                color="white" if value >= 0.55 else "black",
            )
    ax.set_xticks(range(len(CATEGORIES)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(CATEGORIES)), labels)
    ax.set_xlabel("Comparison category")
    ax.set_ylabel("Denominator category")
    ax.set_title("Directional exact-structure overlap among the eleven original positive sets", fontsize=19, weight="bold")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03)
    colorbar.set_label("Fraction of row-category structures also assigned to column category")
    overlap_figure = PAPER / "figures" / "figureS2_directional_category_overlap.png"
    fig.savefig(overlap_figure, dpi=300, facecolor="white")
    plt.close(fig)

    all_structures = set().union(*sets.values())
    multiplicity_counts: dict[int, int] = {}
    for smiles in all_structures:
        multiplicity = sum(smiles in values for values in sets.values())
        multiplicity_counts[multiplicity] = multiplicity_counts.get(multiplicity, 0) + 1
    multiplicity_rows = [
        {
            "category_assignment_count": assignment_count,
            "unique_structure_count": structure_count,
            "fraction_of_all_unique_structures": structure_count / len(all_structures),
        }
        for assignment_count, structure_count in sorted(multiplicity_counts.items())
    ]
    multiplicity_path = PAPER / "supporting_information_category_multiplicity.csv"
    with multiplicity_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(multiplicity_rows[0]))
        writer.writeheader()
        writer.writerows(multiplicity_rows)

    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    bars = ax.bar(
        [row["category_assignment_count"] for row in multiplicity_rows],
        [row["unique_structure_count"] for row in multiplicity_rows],
        color="#4C78A8",
        edgecolor="#333333",
    )
    for bar, row in zip(bars, multiplicity_rows):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(multiplicity_counts.values()) * 0.012,
            f"{int(row['unique_structure_count']):,}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    ax.set_yscale("log")
    ax.set_xticks(sorted(multiplicity_counts))
    ax.set_xlabel("Number of original categories assigned to a structure")
    ax.set_ylabel("Unique structure count (log scale)")
    ax.set_title("Multiplicity of category assignments across unique structures", fontsize=19, weight="bold")
    ax.grid(axis="y", alpha=0.25)
    multiplicity_figure = PAPER / "figures" / "figureS3_category_assignment_multiplicity.png"
    fig.savefig(multiplicity_figure, dpi=300, facecolor="white")
    plt.close(fig)
    return pairwise_path, overlap_figure, multiplicity_path, multiplicity_figure


def add_csv_table(document: Document, path: Path, columns: list[tuple[str, str]]) -> None:
    rows = read_csv(path)
    table = document.add_table(rows=1, cols=len(columns))
    table.style = "Table Grid"
    for cell, (_key, label) in zip(table.rows[0].cells, columns):
        cell.text = label
    for row in rows:
        cells = table.add_row().cells
        for cell, (key, _label) in zip(cells, columns):
            value = row.get(key, "")
            if key.startswith("fraction") or key == "jaccard_similarity":
                value = f"{100 * float(value):.1f}%" if value else ""
            elif "auc_delta" in key:
                value = f"{float(value):+.3f}" if value else ""
            cell.text = value


def build_supporting_document(
    score_selection: Path,
    screening_figure: Path,
    pairwise_path: Path,
    overlap_figure: Path,
    multiplicity_path: Path,
    multiplicity_figure: Path,
    external_database_path: Path,
    combined_test_path: Path,
    combined_test_figure: Path,
) -> Path:
    document = Document()
    normal = document.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(10)
    title = document.add_heading("Supporting Information", 0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    document.add_paragraph(
        "Chemical Category Scoring across Broad Product Classes Using Molecular Descriptors and Structural Patterns"
    ).alignment = WD_ALIGN_PARAGRAPH.CENTER
    document.add_heading("Candidate-score screening", level=1)
    document.add_paragraph(
        "Scoring-function development was attempted for all eleven original categories. Table S1 reports the "
        "exact-overlap-controlled screening result and final disposition of every candidate."
    )
    document.add_paragraph(
        "Table S1. Publication-screening results for all attempted category scoring functions."
    )
    add_csv_table(
        document,
        score_selection,
        [
            ("category", "Category"),
            ("own_threshold_response", "Own response"),
            ("worst_exact_overlap_excluded_source", "Worst foreign source"),
            ("worst_source_threshold_response", "Worst response"),
            ("exact_overlap_excluded_hard_auc", "Hard AUC"),
            ("final_disposition", "Disposition"),
        ],
    )
    document.add_picture(str(screening_figure), width=Inches(6.7))
    document.add_paragraph(
        "Figure S1. Own-category threshold response and worst exact-overlap-excluded foreign-category response "
        "used in publication screening."
    )
    document.add_heading("Exact structural overlap among original categories", level=1)
    document.add_paragraph(
        "Canonical SMILES were deduplicated within each original positive set. For every category pair, Table S2 "
        "reports the exact shared-structure count, the fraction of each category represented by the intersection, "
        "and Jaccard similarity. Figure S2 displays directional overlap, with each cell using the row category as "
        "the denominator."
    )
    document.add_paragraph("Table S2. Pairwise exact-structure overlap among the eleven original positive sets.")
    add_csv_table(
        document,
        pairwise_path,
        [
            ("category_a", "Category A"),
            ("category_b", "Category B"),
            ("shared_structure_count", "Shared"),
            ("fraction_of_category_a", "% of A"),
            ("fraction_of_category_b", "% of B"),
            ("jaccard_similarity", "Jaccard"),
        ],
    )
    document.add_picture(str(overlap_figure), width=Inches(6.7))
    document.add_paragraph(
        "Figure S2. Directional exact-structure overlap among the eleven original positive sets. Values are the "
        "percentage of structures in the row category also assigned to the column category."
    )
    document.add_heading("Multiplicity of category assignments", level=1)
    document.add_paragraph(
        "Assignment multiplicity was counted over the union of unique canonical structures from all eleven positive sets."
    )
    document.add_paragraph("Table S3. Number of original category assignments per unique structure.")
    add_csv_table(
        document,
        multiplicity_path,
        [
            ("category_assignment_count", "Category assignments"),
            ("unique_structure_count", "Unique structures"),
            ("fraction_of_all_unique_structures", "Fraction of union"),
        ],
    )
    document.add_picture(str(multiplicity_figure), width=Inches(6.7))
    document.add_paragraph(
        "Figure S3. Multiplicity of original category assignments across unique structures. The vertical axis uses "
        "a logarithmic scale."
    )
    document.add_heading("External databases and combined-positive-set experiment", level=1)
    document.add_paragraph(
        "Table S4 lists every external database considered for category enrichment or external comparison. Resolved "
        "structures were canonicalized and exact PubChem-category overlaps were counted before defining additions."
    )
    document.add_paragraph("Table S4. External databases considered for positive-set enrichment.")
    add_csv_table(
        document,
        external_database_path,
        [
            ("category", "Category"),
            ("database", "External database"),
            ("citation", "Citation"),
            ("candidate_records", "Candidates"),
            ("resolved_records", "Resolved"),
            ("unique_structures", "Unique"),
            ("new_vs_pubchem", "New vs PubChem"),
            ("use_in_final_test", "Use"),
        ],
    )
    document.add_paragraph(
        "The final experiment added external structures to the PubChem positives within each training fold and rebuilt "
        "the corresponding comparison background after exact target-overlap removal. External additions were divided "
        "into three molecular-hash folds. Table S5 reports candidate-minus-baseline AUC changes on held-out external "
        "positives versus held-out hard cross-category structures and on the original PubChem benchmark."
    )
    document.add_paragraph("Table S5. Three-fold combined-source positive-set rebuilding results.")
    add_csv_table(
        document,
        combined_test_path,
        [
            ("category", "Category"),
            ("pubchem_positive_count", "PubChem positives"),
            ("external_new_count", "External additions"),
            ("mean_external_holdout_auc_delta", "Mean external ΔAUC"),
            ("minimum_external_holdout_auc_delta", "Minimum external ΔAUC"),
            ("mean_original_auc_delta", "Mean original ΔAUC"),
            ("promotion_gate", "Promotion gate"),
        ],
    )
    document.add_picture(str(combined_test_figure), width=Inches(6.7))
    document.add_paragraph(
        "Figure S4. Candidate-minus-baseline AUC changes after combined PubChem and external positive-set rebuilding. "
        "Bars show mean held-out external and original-benchmark changes across three molecular-hash folds."
    )
    output = PAPER / "supporting_information_overlap_analysis.docx"
    document.save(output)
    return output


def build_external_database_assets() -> tuple[Path, Path, Path]:
    rows = [
        {"category": "endocrine_disruptors", "database": "DEDuCT v3 categories I–III", "citation": "[26]", "candidate_records": 704, "resolved_records": 704, "unique_structures": 64, "new_vs_pubchem": 64, "use_in_final_test": "yes"},
        {"category": "flavor_fragrance", "database": "EU Union List of Flavouring Substances", "citation": "[30]", "candidate_records": 2505, "resolved_records": 2264, "unique_structures": 2250, "new_vs_pubchem": 407, "use_in_final_test": "yes"},
        {"category": "pesticides", "database": "Health Canada PMRA PPID", "citation": "[27]", "candidate_records": 545, "resolved_records": 228, "unique_structures": 15, "new_vs_pubchem": 15, "use_in_final_test": "yes"},
        {"category": "surfactants", "database": "US EPA Safer Chemical Ingredients List", "citation": "[31]", "candidate_records": 355, "resolved_records": 181, "unique_structures": 174, "new_vs_pubchem": 162, "use_in_final_test": "yes"},
        {"category": "animal_drugs", "database": "Health Canada Drug Product Database", "citation": "[34]", "candidate_records": 443, "resolved_records": 290, "unique_structures": 284, "new_vs_pubchem": 183, "use_in_final_test": "yes"},
        {"category": "human_drugs", "database": "DrugCentral", "citation": "[33]", "candidate_records": 4099, "resolved_records": 4099, "unique_structures": 1633, "new_vs_pubchem": 1633, "use_in_final_test": "yes"},
        {"category": "human_drugs", "database": "Health Canada Drug Product Database", "citation": "[34]", "candidate_records": 1887, "resolved_records": 14, "unique_structures": 14, "new_vs_pubchem": 13, "use_in_final_test": "partial cached resolution"},
        {"category": "food_additives", "database": "Health Canada Lists of Permitted Food Additives", "citation": "[32]", "candidate_records": 498, "resolved_records": 239, "unique_structures": 236, "new_vs_pubchem": 85, "use_in_final_test": "yes"},
        {"category": "solvents", "database": "US EPA Safer Chemical Ingredients List", "citation": "[31]", "candidate_records": 117, "resolved_records": 101, "unique_structures": 100, "new_vs_pubchem": 75, "use_in_final_test": "yes"},
        {"category": "cosmetics", "database": "California Safe Cosmetics Program", "citation": "[35]", "candidate_records": 113, "resolved_records": 4, "unique_structures": 4, "new_vs_pubchem": 0, "use_in_final_test": "no new structures"},
        {"category": "food_contact_substances", "database": "FSA/FSS regulated-products register", "citation": "[36]", "candidate_records": 763, "resolved_records": 0, "unique_structures": 0, "new_vs_pubchem": 0, "use_in_final_test": "unresolved"},
    ]
    database_path = PAPER / "supporting_information_external_databases.csv"
    with database_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = json.loads(
        (ROOT / "results/combined_external_positive_rebuild/summary.json").read_text(encoding="utf-8")
    )
    test_rows: list[dict[str, object]] = []
    external_counts = {
        row["category"]: int(row["new_vs_pubchem"])
        for row in rows
        if row["use_in_final_test"] == "yes"
    }
    external_counts["human_drugs"] = 1639
    for category, values in summary.items():
        folds = values["folds"]
        test_rows.append(
            {
                "category": category,
                "pubchem_positive_count": folds[0]["pubchem_positive_count"],
                "external_new_count": external_counts[category],
                "mean_external_holdout_auc_delta": values["mean_external_holdout_auc_delta"],
                "minimum_external_holdout_auc_delta": values["minimum_external_holdout_auc_delta"],
                "mean_original_auc_delta": values["mean_original_auc_delta"],
                "promotion_gate": str(values["promotion_gate"]).lower(),
            }
        )
    test_path = PAPER / "supporting_information_combined_positive_rebuild.csv"
    with test_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(test_rows[0]))
        writer.writeheader()
        writer.writerows(test_rows)

    x = np.arange(len(test_rows))
    width = 0.36
    fig, ax = plt.subplots(figsize=(14, 8), constrained_layout=True)
    ax.bar(
        x - width / 2,
        [float(row["mean_external_holdout_auc_delta"]) for row in test_rows],
        width,
        label="Held-out external ΔAUC",
        color="#4C78A8",
    )
    ax.bar(
        x + width / 2,
        [float(row["mean_original_auc_delta"]) for row in test_rows],
        width,
        label="Original-benchmark ΔAUC",
        color="#F58518",
    )
    ax.axhline(0, color="black", linewidth=1)
    ax.set_xticks(x, [DISPLAY.get(row["category"], row["category"]).replace("_", " ") for row in test_rows], rotation=40, ha="right")
    ax.set_ylabel("Candidate minus baseline AUC")
    ax.set_title("Combined PubChem and external positive-set rebuilding", fontsize=20, weight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    figure_path = PAPER / "figures/figureS4_combined_positive_rebuild.png"
    fig.savefig(figure_path, dpi=300, facecolor="white")
    plt.close(fig)
    return database_path, test_path, figure_path


def main() -> None:
    rows = read_csv(ANALYSIS / "publication_specificity_revalidation.csv")
    animal_human = [
        {
            "category": "animal_drugs",
            "own_threshold_response": "",
            "worst_source_threshold_response": "",
            "exact_overlap_excluded_hard_auc": "",
            "publication_decision": "exclude: QED overlap",
        },
        {
            "category": "human_drugs",
            "own_threshold_response": "",
            "worst_source_threshold_response": "",
            "exact_overlap_excluded_hard_auc": "",
            "publication_decision": "exclude: QED overlap",
        },
    ]
    output_rows = animal_human + rows
    table_path = PAPER / "supporting_information_score_selection.csv"
    with table_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = [
            "category",
            "own_threshold_response",
            "worst_exact_overlap_excluded_source",
            "worst_source_threshold_response",
            "exact_overlap_excluded_hard_auc",
            "publication_decision",
            "final_disposition",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in output_rows:
            category = row["category"]
            decision = row["publication_decision"]
            if category in {"flavoring_agents", "fragrances"}:
                disposition = "merge as flavor and fragrance"
            elif decision == "retain":
                disposition = "retain"
            else:
                disposition = "exclude"
            writer.writerow({field: row.get(field, "") for field in fields[:-1]} | {"final_disposition": disposition})

    plot_rows = [row for row in rows if row["category"] not in {"flavoring_agents", "fragrances"}]
    categories = [row["category"] for row in plot_rows]
    own = np.asarray([float(row["own_threshold_response"]) for row in plot_rows])
    worst = np.asarray([float(row["worst_source_threshold_response"]) for row in plot_rows])
    x = np.arange(len(categories)); width = 0.36
    fig, ax = plt.subplots(figsize=(14, 8), constrained_layout=True)
    ax.bar(x - width / 2, own, width, color="#2C7FB8", label="Own positive set")
    ax.bar(x + width / 2, worst, width, color="#F28E2B", label="Worst exact-overlap-excluded source")
    for index, row in enumerate(plot_rows):
        retained = row["publication_decision"] == "retain"
        ax.text(index, max(own[index], worst[index]) + 0.035, "Retain" if retained else "Exclude", ha="center", weight="bold", color="#1B7837" if retained else "#B2182B")
    ax.set_title("Publication screening of attempted category scores", fontsize=21, weight="bold")
    ax.set_ylabel("Fraction at or above frozen threshold", fontsize=15)
    ax.set_xticks(x, [DISPLAY[c] for c in categories], rotation=40, ha="right")
    ax.set_ylim(0, 1.12); ax.grid(axis="y", alpha=0.25); ax.legend(fontsize=12)
    fig.savefig(PAPER / "figures" / "figureS1_publication_score_screening.png", dpi=300, facecolor="white")
    plt.close(fig)
    pairwise_path, overlap_figure, multiplicity_path, multiplicity_figure = build_overlap_assets()
    external_database_path, combined_test_path, combined_test_figure = build_external_database_assets()
    supporting_document = build_supporting_document(
        table_path,
        PAPER / "figures" / "figureS1_publication_score_screening.png",
        pairwise_path,
        overlap_figure,
        multiplicity_path,
        multiplicity_figure,
        external_database_path,
        combined_test_path,
        combined_test_figure,
    )
    print(table_path)
    print(PAPER / "figures" / "figureS1_publication_score_screening.png")
    print(pairwise_path)
    print(overlap_figure)
    print(multiplicity_path)
    print(multiplicity_figure)
    print(supporting_document)
    print(external_database_path)
    print(combined_test_path)
    print(combined_test_figure)


if __name__ == "__main__":
    main()
