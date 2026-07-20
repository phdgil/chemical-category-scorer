from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
RESULTS_DIR = ROOT_DIR / 'results'
PAPER_DIR = ROOT_DIR / 'paper'

FINAL_CSV = RESULTS_DIR / 'final_scoring_performance_table.csv'
FINAL_JSON = RESULTS_DIR / 'final_scoring_performance_table.json'
PAPER_FINAL_CSV = PAPER_DIR / 'final_scoring_performance_table.csv'
QED_CSV = RESULTS_DIR / 'qed_baseline_comparison' / 'qed_baseline_comparison.csv'

QED_FIELDS = [
    'qed_baseline_auc',
    'qed_baseline_balanced_accuracy',
    'delta_auc_vs_qed',
    'delta_balanced_accuracy_vs_qed',
    'beats_qed_auc',
    'beats_qed_balanced_accuracy',
    'reportable_vs_qed',
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def enrich_rows(final_rows: list[dict[str, str]], qed_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    qed_by_id = {row['output_id']: row for row in qed_rows}
    enriched: list[dict[str, str]] = []
    for row in final_rows:
        merged = dict(row)
        if row.get('output_level') == 'category' and row.get('output_id') in qed_by_id:
            qed = qed_by_id[row['output_id']]
            merged.update(
                {
                    'qed_baseline_auc': qed.get('qed_auc', ''),
                    'qed_baseline_balanced_accuracy': qed.get('qed_balanced_accuracy', ''),
                    'delta_auc_vs_qed': qed.get('delta_auc_vs_qed', ''),
                    'delta_balanced_accuracy_vs_qed': qed.get('delta_balanced_accuracy_vs_qed', ''),
                    'beats_qed_auc': qed.get('beats_qed_auc', ''),
                    'beats_qed_balanced_accuracy': qed.get('beats_qed_balanced_accuracy', ''),
                    'reportable_vs_qed': qed.get('reportable_vs_qed', ''),
                }
            )
        else:
            for field in QED_FIELDS:
                merged.setdefault(field, '')
                merged[field] = ''
        enriched.append(merged)
    return enriched


def main() -> None:
    final_rows = read_csv(FINAL_CSV)
    qed_rows = read_csv(QED_CSV)
    enriched = enrich_rows(final_rows, qed_rows)
    write_csv(FINAL_CSV, enriched)
    write_csv(PAPER_FINAL_CSV, enriched)
    FINAL_JSON.write_text(json.dumps(enriched, indent=2, ensure_ascii=False), encoding='utf-8')
    print(f'Enriched {len(enriched)} final rows with QED baselines')


if __name__ == '__main__':
    main()
