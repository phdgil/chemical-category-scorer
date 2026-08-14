from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median

from app.algorithm_score_engine import score_smiles
from app.audit_external_validation_overlap import run_audit


DATASETS = {
    "animal_drugs": "animal_drugs_canada_dpd.csv",
    "cosmetics": "cosmetics_california_current_reports.csv",
    "flavoring_agents": "flavoring_agents_eu.csv",
    "food_additives": "food_additives_canada.csv",
    "food_contact_substances": "food_contact_gb_authorised.csv",
    "human_drugs": "human_drugs_drugcentral.csv",
    "pesticides": "pesticides_canada_pmra.csv",
}

MODEL_IDS = {
    "animal_drugs": "final_animal_drugs",
    "cosmetics": "final_cosmetics",
    "endocrine_disruptors": "han_endocrine_disruptors",
    "flavoring_agents": "final_flavoring_agents",
    "food_additives": "final_food_additives",
    "food_contact_substances": "final_food_contact_substances",
    "human_drugs": "final_human_drugs",
    "pesticides": "final_pesticides",
}

SOURCE_LABELS = {
    "animal_drugs": "Health Canada DPD (Veterinary)",
    "cosmetics": "California Safe Cosmetics Program",
    "endocrine_disruptors": "DEDuCT v3 I-III",
    "flavoring_agents": "EU Union List of Flavouring Substances",
    "food_additives": "Health Canada permitted-additive lists",
    "food_contact_substances": "FSA/FSS GB food-contact register",
    "human_drugs": "DrugCentral",
    "pesticides": "Health Canada PMRA PPID",
}

SOURCE_INDEPENDENCE = {
    "animal_drugs": "independent regulator; annotation-lineage independence not proven",
    "cosmetics": "known PubChem source overlap; external-structure audit only",
    "flavoring_agents": "known PubChem source-family overlap; external-structure audit only",
    "food_additives": "independent regulator; annotation-lineage independence not proven",
    "food_contact_substances": "plausibly independent regulator; annotation-lineage independence not proven",
    "human_drugs": "independently curated database; annotation-lineage independence not proven",
    "pesticides": "independent regulator; annotation-lineage independence not proven",
}

EVIDENCE_SCOPE = {
    "animal_drugs": "external positive comparison",
    "cosmetics": "source-overlapping consistency comparison",
    "endocrine_disruptors": "external positive comparison",
    "flavoring_agents": "source-overlapping consistency comparison",
    "food_additives": "external positive comparison",
    "food_contact_substances": "external positive comparison",
    "human_drugs": "external positive comparison",
    "pesticides": "external positive comparison",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def score_candidates(category: str, rows: list[dict[str, str]]) -> dict[str, object]:
    scored: list[float] = []
    recovered = 0
    invalid = 0
    threshold: float | None = None
    for row in rows:
        result = score_smiles(row.get("standardized_smiles") or row.get("SMILES", ""), model_id=MODEL_IDS[category])
        if not result.valid:
            invalid += 1
            continue
        scored.append(float(result.score))
        threshold = float(result.threshold)
        recovered += int(result.score >= result.threshold)
    return {
        "scored_true_external": len(scored),
        "invalid_at_scoring": invalid,
        "frozen_threshold": threshold,
        "recovered_above_threshold": recovered,
        "positive_recovery_fraction": recovered / len(scored) if scored else None,
        "median_score": median(scored) if scored else None,
        "minimum_score": min(scored) if scored else None,
        "maximum_score": max(scored) if scored else None,
    }


def run(resolved_dir: Path, inputs_dir: Path, output_dir: Path) -> list[dict[str, object]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []

    endocrine_path = output_dir.parent / "deduct_v3_endocrine" / "endocrine_disruptors_true_external_candidates_deduct_I-III.csv"
    if endocrine_path.is_file():
        endocrine_rows = read_csv(endocrine_path)
        row: dict[str, object] = {
            "category": "endocrine_disruptors",
            "external_source": SOURCE_LABELS["endocrine_disruptors"],
            "evidence_scope": EVIDENCE_SCOPE["endocrine_disruptors"],
            "raw_candidates": 704,
            "resolved_structures": 704,
            "overlap_excluded": 639,
            "duplicate_external_excluded": 1,
            "unresolved": 0,
            "true_external_candidates": len(endocrine_rows),
            "source_independence": "literature-curated labels; molecule overlap removed",
        }
        row.update(score_candidates("endocrine_disruptors", endocrine_rows))
        summaries.append(row)

    for category, filename in DATASETS.items():
        resolved_path = resolved_dir / filename
        if not resolved_path.is_file():
            continue
        resolution_rows = read_csv(resolved_path)
        resolved_count = sum(row.get("resolution_status") == "resolved" for row in resolution_rows)
        audit_dir = output_dir / category
        audit_summary = run_audit(
            external_file=resolved_path,
            inputs_dir=inputs_dir,
            target_category=category,
            smiles_column="SMILES",
            identifier_column="external_id",
            output_dir=audit_dir,
        )
        status = audit_summary["status_counts"]
        true_rows = read_csv(audit_dir / f"{category}_true_external_candidates.csv")
        row = {
            "category": category,
            "external_source": SOURCE_LABELS[category],
            "evidence_scope": EVIDENCE_SCOPE[category],
            "raw_candidates": len(resolution_rows),
            "resolved_structures": resolved_count,
            "overlap_excluded": status.get("overlap_excluded", 0),
            "duplicate_external_excluded": status.get("duplicate_external_excluded", 0),
            "unresolved": status.get("unresolved", 0),
            "true_external_candidates": status.get("true_external_candidate", 0),
            "source_independence": SOURCE_INDEPENDENCE[category],
        }
        row.update(score_candidates(category, true_rows))
        summaries.append(row)

    summaries.sort(key=lambda row: str(row["category"]))
    write_csv(output_dir / "external_validation_summary.csv", summaries)
    (output_dir / "external_validation_summary.json").write_text(
        json.dumps(summaries, indent=2), encoding="utf-8"
    )
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolved-dir", type=Path, required=True)
    parser.add_argument("--inputs-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.resolved_dir, args.inputs_dir, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
