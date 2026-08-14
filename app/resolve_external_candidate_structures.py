from __future__ import annotations

import argparse
import csv
import json
import subprocess
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def load_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_cache(path: Path, cache: dict[str, dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(path)


def request_json(url: str, *, attempts: int = 6) -> dict[str, Any]:
    delay = 1.0
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": "chemical-category-external-validation/1.0"})
            with urlopen(request, timeout=60) as response:
                result = json.loads(response.read().decode("utf-8"))
            time.sleep(0.25)
            return result
        except HTTPError as exc:
            if exc.code == 404:
                return {"not_found": True}
            if exc.code not in {429, 500, 502, 503, 504} or attempt + 1 == attempts:
                return {"error": f"HTTP {exc.code}"}
        except (TimeoutError, URLError, json.JSONDecodeError) as exc:
            if attempt + 1 == attempts:
                return {"error": f"{type(exc).__name__}: {exc}"}
        time.sleep(delay)
        delay = min(delay * 2, 30)
    return {"error": "request attempts exhausted"}


def resolve_with_nci(query: str) -> dict[str, Any]:
    """Use the NCI Chemical Identifier Resolver through curl."""
    encoded = quote(query, safe="")
    url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/smiles"
    completed = subprocess.run(
        [
            "curl.exe",
            "--silent",
            "--show-error",
            "--retry",
            "1",
            "--retry-delay",
            "1",
            "--retry-all-errors",
            "--retry-max-time",
            "30",
            "--connect-timeout",
            "10",
            "--max-time",
            "20",
            "--write-out",
            "\\n%{http_code}",
            url,
        ],
        capture_output=True,
        text=True,
        timeout=45,
        check=False,
    )
    lines = completed.stdout.rstrip().splitlines()
    status_code = lines[-1] if lines else ""
    smiles = "\\n".join(lines[:-1]).strip() if len(lines) > 1 else ""
    if completed.returncode != 0:
        return {
            "status": "request_error",
            "query": query,
            "error": completed.stderr.strip() or f"curl exit {completed.returncode}",
        }
    if status_code == "404":
        return {"status": "not_found", "query": query, "error": "HTTP 404"}
    if status_code != "200":
        return {"status": "request_error", "query": query, "error": f"HTTP {status_code}"}
    if not smiles or "resolver" in smiles.lower():
        return {"status": "not_found", "query": query}
    return {
        "status": "resolved",
        "query": query,
        "smiles": smiles,
        "resolution_source": "NCI_CIR",
    }


def resolve_query(query: str, *, require_unique: bool) -> dict[str, Any]:
    encoded = quote(query, safe="")
    if require_unique:
        cid_result = request_json(f"{PUG_BASE}/compound/name/{encoded}/cids/JSON")
        if cid_result.get("not_found"):
            return {"status": "not_found", "query": query}
        if cid_result.get("error"):
            return {"status": "request_error", "query": query, "error": cid_result["error"]}
        cids = cid_result.get("IdentifierList", {}).get("CID", [])
        if len(cids) != 1:
            return {"status": "ambiguous", "query": query, "candidate_cids": cids}
        cid = str(cids[0])
        property_result = request_json(
            f"{PUG_BASE}/compound/cid/{cid}/property/SMILES,ConnectivitySMILES,InChIKey,IUPACName/JSON"
        )
    else:
        property_result = request_json(
            f"{PUG_BASE}/compound/name/{encoded}/property/SMILES,ConnectivitySMILES,InChIKey,IUPACName/JSON"
        )

    if property_result.get("not_found"):
        return {"status": "not_found", "query": query}
    if property_result.get("error"):
        return {"status": "request_error", "query": query, "error": property_result["error"]}
    properties = property_result.get("PropertyTable", {}).get("Properties", [])
    if len(properties) != 1:
        return {
            "status": "ambiguous" if properties else "not_found",
            "query": query,
            "candidate_count": len(properties),
        }
    prop = properties[0]
    smiles = prop.get("SMILES") or prop.get("CanonicalSMILES") or prop.get("ConnectivitySMILES")
    if not smiles or not prop.get("InChIKey"):
        return {"status": "unresolved_structure", "query": query, "cid": prop.get("CID")}
    return {
        "status": "resolved",
        "query": query,
        "cid": prop.get("CID"),
        "smiles": smiles,
        "connectivity_smiles": prop.get("ConnectivitySMILES", ""),
        "inchikey": prop.get("InChIKey", ""),
        "iupac_name": prop.get("IUPACName", ""),
        "resolution_source": "PubChem_PUG_REST",
    }


def resolve_candidate_file(
    input_path: Path,
    output_path: Path,
    cache_path: Path,
    *,
    checkpoint_every: int = 25,
    prefer_nci: bool = False,
    cache_only: bool = False,
) -> dict[str, int]:
    cache = load_cache(cache_path)
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    counts: dict[str, int] = {}
    output_rows: list[dict[str, Any]] = []

    for index, row in enumerate(rows, start=1):
        cas = (row.get("cas") or "").strip()
        name = (row.get("external_name") or "").strip()
        query = cas or name
        method = "cas" if cas else "name"
        cache_key = f"{method}:{query}"
        if not query:
            resolution = {"status": "missing_identifier", "query": ""}
        elif cache_only:
            resolution = cache.get(
                cache_key,
                {"status": "not_attempted", "query": query, "error": "No cached resolution"},
            )
        elif cache_key in cache and cache[cache_key].get("status") not in {"request_error", "not_found"}:
            resolution = cache[cache_key]
        else:
            if prefer_nci:
                resolution = resolve_with_nci(query)
            else:
                resolution = resolve_query(query, require_unique=(method == "name"))
                if resolution.get("status") in {"request_error", "not_found"}:
                    resolution = resolve_with_nci(query)
            cache[cache_key] = resolution
        status = str(resolution.get("status", "unknown"))
        counts[status] = counts.get(status, 0) + 1
        output = dict(row)
        output.update(
            {
                "resolution_query": query,
                "resolution_method": method,
                "resolution_status": status,
                "pubchem_cid": resolution.get("cid", ""),
                "SMILES": resolution.get("smiles", ""),
                "connectivity_smiles": resolution.get("connectivity_smiles", ""),
                "resolved_inchikey": resolution.get("inchikey", ""),
                "resolved_iupac_name": resolution.get("iupac_name", ""),
                "resolution_error": resolution.get("error", ""),
                "resolution_source": resolution.get("resolution_source", ""),
                "candidate_cids": ";".join(map(str, resolution.get("candidate_cids", []))),
            }
        )
        output_rows.append(output)
        if index % checkpoint_every == 0:
            save_cache(cache_path, cache)
            print(f"{input_path.name}: {index}/{len(rows)} resolved={counts.get('resolved', 0)}", flush=True)

    save_cache(cache_path, cache)
    fields: list[str] = []
    seen: set[str] = set()
    for row in output_rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"completed {input_path.name}: {json.dumps(counts, sort_keys=True)}", flush=True)
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resolve external candidate identifiers to structures using PubChem PUG REST.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--prefer-nci", action="store_true")
    parser.add_argument("--cache-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    resolve_candidate_file(
        args.input,
        args.output,
        args.cache,
        prefer_nci=args.prefer_nci,
        cache_only=args.cache_only,
    )


if __name__ == "__main__":
    main()
