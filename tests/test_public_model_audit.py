from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from app import public_model_audit as audit


class PublicModelAuditTests(unittest.TestCase):
    def test_public_model_roles_are_release_ordered_and_public_only(self) -> None:
        model_ids = audit.public_model_ids()

        self.assertEqual(4, len(model_ids))
        self.assertNotIn("final_endocrine_disruptors", model_ids)
        self.assertEqual(["han_endocrine_disruptors"], [model_id for model_id in model_ids if audit.get_model_role(model_id) == audit.AUXILIARY_HAZARD_ROLE])
        self.assertEqual(3, len(audit.product_model_ids()))

    def test_probe_audit_matches_committed_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            panel = tmp_path / "panel.csv"
            scores = tmp_path / "scores.csv"
            summary = tmp_path / "summary.csv"

            audit.write_default_probe_panel(panel)
            audit.run_audit(panel, scores, summary)

            self.assertEqual(audit.DEFAULT_PROBE_PANEL.read_text(encoding="utf-8"), panel.read_text(encoding="utf-8"))
            self.assertEqual(audit.DEFAULT_AUDIT_SCORES.read_text(encoding="utf-8"), scores.read_text(encoding="utf-8"))
            self.assertEqual(audit.DEFAULT_AUDIT_SUMMARY.read_text(encoding="utf-8"), summary.read_text(encoding="utf-8"))

    def test_audit_summary_separates_product_and_auxiliary_signal(self) -> None:
        scores, summaries = audit.build_audit_rows(audit.DEFAULT_PROBE_PANEL)
        by_probe = {row["probe_id"]: row for row in summaries}
        score_fields = set(scores[0])

        ddt = by_probe["ddt"]
        self.assertIn("release_model_order", score_fields)
        self.assertIn("score_rank_within_role", score_fields)
        self.assertNotIn("model_rank", score_fields)
        self.assertEqual("han_endocrine_disruptors", ddt["auxiliary_endocrine_model_id"])
        self.assertNotEqual("han_endocrine_disruptors", ddt["top_product_model_id"])
        self.assertIn("pesticides", ddt["top_three_product_categories"])
        self.assertTrue(ddt["likely_product_count"].isdigit())

    def test_pattern_overlap_has_shared_fixed_smarts_and_docs_limitation(self) -> None:
        rows = audit.build_pattern_overlap_rows()
        shared = [row for row in rows if row["shared_across_models"] == "true"]
        shared_names = " ".join(row["pattern_names"] for row in shared)

        self.assertEqual(0, len(shared))
        self.assertEqual("", shared_names)
        markdown = audit.DEFAULT_PATTERN_AUDIT.read_text(encoding="utf-8")
        self.assertIn("not a labeled classification benchmark", markdown)
        self.assertIn("makes no accuracy or BRICS-improvement claim", markdown)
        self.assertIn("## Pattern-unit panel data quality", markdown)
        self.assertIn("Rows: 65", markdown)
        self.assertIn("Columns: 12", markdown)
        self.assertIn("Unique test SMILES: 25", markdown)
        self.assertIn("Repeated probe SMILES are expected", markdown)

    def test_app_test_smiles_quality_is_computed(self) -> None:
        quality = audit.app_test_smiles_quality()

        self.assertEqual(65, quality["row_count"])
        self.assertEqual(12, quality["column_count"])
        self.assertEqual(25, quality["unique_test_smiles_count"])
        self.assertEqual(0, quality["missing_cell_count"])
        self.assertEqual(0, quality["exact_duplicate_row_count"])
        self.assertEqual(0, quality["invalid_test_smiles_count"])
        self.assertEqual(0, quality["invalid_pattern_smarts_count"])

    def test_candidate_extraction_is_deterministic_and_handles_invalid_smiles(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            input_csv = tmp_path / "input.csv"
            first = tmp_path / "first.csv"
            second = tmp_path / "second.csv"
            summary = tmp_path / "summary.csv"
            input_csv.write_text(
                "category,SMILES\n"
                "pesticides,Clc1ccccc1\n"
                "pesticides,CC(=O)OC\n"
                "flavoring_agents,COc1cc(C=O)ccc1O\n"
                "flavoring_agents,CC(=O)OC\n"
                "cosmetics,COC(=O)/C=C/c1ccccc1\n"
                "cosmetics,CC(=O)OC\n"
                "bad,not-a-smiles\n",
                encoding="utf-8",
            )

            audit.write_pattern_candidates(first, summary_csv=summary, input_csv=input_csv, min_prevalence=1, min_enrichment=1.0)
            audit.write_pattern_candidates(second, input_csv=input_csv, min_prevalence=1, min_enrichment=1.0)

            self.assertEqual(first.read_text(encoding="utf-8"), second.read_text(encoding="utf-8"))
            with first.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            with summary.open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            self.assertTrue(rows)
            self.assertEqual({"fixed_smarts", "murcko", "brics", "hybrid"}, {row["method"] for row in summary_rows})
            self.assertTrue(all(row["claim_scope"] == "diversity_comparison_only_no_performance_claim" for row in summary_rows))
            self.assertTrue(all("mean_pairwise_jaccard" in row for row in summary_rows))
            self.assertTrue(all("shared_three_plus_unit_fraction" in row for row in summary_rows))
            self.assertTrue(all("mean_category_breadth" in row for row in summary_rows))
            fixed_summary = next(row for row in summary_rows if row["method"] == "fixed_smarts")
            self.assertGreater(float(fixed_summary["mean_pairwise_jaccard"]), 0.0)
            self.assertGreater(float(fixed_summary["shared_three_plus_unit_fraction"]), 0.0)
            self.assertGreaterEqual(float(fixed_summary["mean_category_breadth"]), 1.0)
            self.assertTrue(all(row["invalid_smiles_skipped"] == "1" for row in rows))
            brics_units = {row["unit"] for row in rows if row["method"] == "brics"}
            self.assertIn("O=Cc1ccc(O)cc1", brics_units)
            self.assertFalse(any("*" in unit for unit in brics_units))

    def test_fragment_helpers_cover_murcko_brics_canonicalization(self) -> None:
        aspirin = audit.canonical_smiles("CC(=O)Oc1ccccc1C(=O)O")

        self.assertEqual("CC(=O)Oc1ccccc1C(=O)O", aspirin)
        self.assertEqual("", audit.canonical_smiles("not-a-smiles"))
        self.assertIn("c1ccccc1", audit.extract_candidate_units(aspirin, "murcko"))
        self.assertTrue(audit.extract_candidate_units(aspirin, "brics"))
        self.assertTrue(audit.extract_candidate_units(aspirin, "hybrid"))

    def test_hash_sampling_is_not_file_order_first_n(self) -> None:
        ordered = ["CCO", "CCCC", "c1ccccc1", "CCN", "CCCl"]
        sampled = audit._hash_sample(ordered, 2)

        self.assertEqual(sampled, audit._hash_sample(list(reversed(ordered)), 2))
        self.assertNotEqual(ordered[:2], sampled)

    def test_missing_final_rebuild_inputs_fail_fast(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            missing_dir = Path(tmp) / "missing-final-inputs"
            with patch.object(audit, "FINAL_REBUILD_INPUTS_DIR", missing_dir):
                with self.assertRaisesRegex(FileNotFoundError, "Missing final rebuild input files"):
                    audit.read_final_rebuild_inputs(limit_per_category=1)


if __name__ == "__main__":
    unittest.main()
