from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

APP_DIR = Path(__file__).resolve().parent
RESULTS_DIR = APP_DIR.parent / "results"
CATALOG_CSV = RESULTS_DIR / "pubchem_hid72_catalog.csv"
OUT_JSON = RESULTS_DIR / "pubchem_named_subtype_sources.json"
OUT_CSV = RESULTS_DIR / "pubchem_named_subtype_sources.csv"

SOURCE_MAP: dict[str, list[dict[str, Any]]] = {
    "cosmetics": [
        {
            "pubchem_source_name": "California Safe Cosmetics Program (CSCP) Product Database",
            "expected_role": "official named cosmetics-related PubChem source for future subtype labeling",
        },
        {
            "pubchem_source_name": "Consumer Product Information Database Classification",
            "expected_role": "product-use classification source with consumer-product terms",
        },
        {
            "pubchem_source_name": "EPA CPDat Classification",
            "expected_role": "consumer-product and exposure-use taxonomy candidate",
        },
    ],
    "food_contact_substances": [
        {
            "pubchem_source_name": "FDA Food Contact Substances (FCS)",
            "expected_role": "official named food-contact source for subtype replacement",
        },
        {
            "pubchem_source_name": "FDA Packaging & Food Contact Substances (FCS)",
            "expected_role": "official packaging/food-contact source mirror",
        },
        {
            "pubchem_source_name": "Food Ontology (FoodOn)",
            "expected_role": "food-related controlled vocabulary for finer naming where applicable",
        },
    ],
    "human_drugs": [
        {
            "pubchem_source_name": "FDA Drug Type and Pharmacologic Classification",
            "expected_role": "official named human-drug subtype source",
            "source_hid": 116,
        },
        {
            "pubchem_source_name": "FDA Pharm Classes",
            "expected_role": "official pharmacologic class source",
            "source_hid": 78,
        },
        {
            "pubchem_source_name": "KEGG: Drug",
            "expected_role": "named drug taxonomy candidate",
        },
    ],
    "animal_drugs": [
        {
            "pubchem_source_name": "ATCvet Classification",
            "expected_role": "official named animal-drug subtype source",
            "source_hid": 136,
        },
        {
            "pubchem_source_name": "FDA Pharm Classes",
            "expected_role": "shared pharmacologic class source where veterinary mapping exists",
            "source_hid": 78,
        },
        {
            "pubchem_source_name": "KEGG: Drug",
            "expected_role": "fallback named drug taxonomy candidate",
        },
    ],
}


def read_catalog() -> list[dict[str, Any]]:
    with CATALOG_CSV.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def build_rows() -> list[dict[str, Any]]:
    catalog = read_catalog()
    by_name = {row["name"]: row for row in catalog}
    rows: list[dict[str, Any]] = []
    for panel_id, specs in SOURCE_MAP.items():
        for spec in specs:
            row = by_name.get(spec["pubchem_source_name"])
            rows.append(
                {
                    "panel_id": panel_id,
                    "pubchem_source_name": spec["pubchem_source_name"],
                    "expected_role": spec["expected_role"],
                    "found_in_catalog": bool(row),
                    "source_hid": spec.get("source_hid", ""),
                    "hnid": row.get("hnid", "") if row else "",
                    "compound_count": row.get("compound_count", "") if row else "",
                    "path": row.get("path", "") if row else "",
                }
            )
    return rows


def main() -> None:
    rows = build_rows()
    OUT_JSON.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["panel_id", "pubchem_source_name", "expected_role", "found_in_catalog", "source_hid", "hnid", "compound_count", "path"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(OUT_JSON)
    print(OUT_CSV)


if __name__ == "__main__":
    main()
