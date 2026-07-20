from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

APP_DIR = Path(__file__).resolve().parent
RESULTS_DIR = APP_DIR.parent / "results"
PANELS_DIR = APP_DIR / "data" / "evidence_panels"
SOURCE_MAP_JSON = RESULTS_DIR / "pubchem_named_subtype_sources.json"
OFFICIAL_DRUG_CATALOG_JSON = RESULTS_DIR / "official_drug_subtype_validation" / "official_drug_subtype_catalog.json"
OUT_JSON = RESULTS_DIR / "pubchem_named_subtype_taxonomy.json"
OUT_CSV = RESULTS_DIR / "pubchem_named_subtype_taxonomy.csv"
PROTOTYPE_PANELS = ["cosmetics", "food_contact_substances"]
OFFICIAL_DRUG_PANELS = ["human_drugs", "animal_drugs"]


def humanize(text: str) -> str:
    return text.replace("_", " ").strip().title()


def load_source_map() -> dict[str, list[dict[str, Any]]]:
    rows = json.loads(SOURCE_MAP_JSON.read_text(encoding="utf-8"))
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row["panel_id"]), []).append(row)
    return grouped


def build_prototype_rows(source_map: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for panel_id in PROTOTYPE_PANELS:
        panel_path = PANELS_DIR / f"{panel_id}.json"
        panel = json.loads(panel_path.read_text(encoding="utf-8"))
        source_candidates = source_map.get(panel_id, [])
        source_names = [row["pubchem_source_name"] for row in source_candidates]
        source_ids = [str(row.get("hnid") or row.get("source_hid") or "") for row in source_candidates if row.get("hnid") or row.get("source_hid")]
        for family in panel.get("families", []):
            family_id = str(family["family_id"])
            family_name = str(family.get("family_name", family_id))
            rows.append(
                {
                    "panel_id": panel_id,
                    "family_id": family_id,
                    "display_name": f"{humanize(panel_id)} Prototype / {humanize(family_name)}",
                    "prototype_family_name": family_name,
                    "member_count": int(family.get("member_count", 0)),
                    "pubchem_official_named_subtype_available": False,
                    "naming_mode": "prototype_family_plus_pubchem_source_reference",
                    "pubchem_source_candidates": "; ".join(source_names),
                    "pubchem_source_hnids": "; ".join(source_ids),
                    "representative_smiles": str(family.get("representative_smiles", "")),
                }
            )
    return rows


def build_official_drug_rows() -> list[dict[str, Any]]:
    if not OFFICIAL_DRUG_CATALOG_JSON.exists():
        return []
    catalog_rows = json.loads(OFFICIAL_DRUG_CATALOG_JSON.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for row in catalog_rows:
        if str(row.get("panel_id", "")) not in OFFICIAL_DRUG_PANELS:
            continue
        if not bool(row.get("selected_for_validation")):
            continue
        subtype_name = str(row.get("subtype_name", ""))
        rows.append(
            {
                "panel_id": str(row["panel_id"]),
                "family_id": f"official_{row['subtype_hnid']}",
                "display_name": f"{humanize(str(row['panel_id']))} Official / {subtype_name}",
                "prototype_family_name": subtype_name,
                "member_count": int(row.get("hid72_overlap_unique_smiles_count", 0)),
                "pubchem_official_named_subtype_available": True,
                "naming_mode": "official_pubchem_subtype",
                "pubchem_source_candidates": str(row.get("source_name", "")),
                "pubchem_source_hnids": f"source_hid={row.get('source_hid', '')}; subtype_hnid={row.get('subtype_hnid', '')}",
                "representative_smiles": str(row.get("representative_smiles", "")),
            }
        )
    return rows


def build_rows() -> list[dict[str, Any]]:
    source_map = load_source_map()
    rows = build_prototype_rows(source_map)
    rows.extend(build_official_drug_rows())
    rows.sort(key=lambda row: (row["panel_id"], row["display_name"]))
    return rows


def main() -> None:
    rows = build_rows()
    OUT_JSON.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "panel_id",
                "family_id",
                "display_name",
                "prototype_family_name",
                "member_count",
                "pubchem_official_named_subtype_available",
                "naming_mode",
                "pubchem_source_candidates",
                "pubchem_source_hnids",
                "representative_smiles",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(OUT_JSON)
    print(OUT_CSV)


if __name__ == "__main__":
    main()
