from __future__ import annotations

import argparse
import base64
import csv
import tempfile
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

if __package__:
    from .algorithm_score_engine import (
        ALL_MODELS_SENTINEL,
        AUXILIARY_HAZARD_ROLE,
        MODELS_DIR,
        PRODUCT_USE_ROLE,
        choose_best_product_result,
        choose_best_result,
        choose_representative_product_result,
        cross_category_specificity,
        get_model_role,
        list_models,
        refresh_model_registry,
        render_molecule_png,
        render_pattern_match_png,
        score_csv,
        score_smiles,
        score_smiles_all,
    )
else:  # pragma: no cover - for direct script execution
    from algorithm_score_engine import (
        ALL_MODELS_SENTINEL,
        AUXILIARY_HAZARD_ROLE,
        MODELS_DIR,
        PRODUCT_USE_ROLE,
        choose_best_product_result,
        choose_best_result,
        choose_representative_product_result,
        cross_category_specificity,
        get_model_role,
        list_models,
        refresh_model_registry,
        render_molecule_png,
        render_pattern_match_png,
        score_csv,
        score_smiles,
        score_smiles_all,
    )

APP_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = APP_DIR / "output"
LOGO_PATH = APP_DIR / "data" / "AAA_logo.png"
DEFAULT_MODEL_ID = "final_pesticides"


def _clean_label(label: str) -> str:
    return str(label).replace(" / Final rebuilt scorer", "")


def _clean_decision(text: str) -> str:
    text = str(text or "").replace("_", " ").strip()
    return text[:1].upper() + text[1:] if text else text


def _humanize_pattern(name: str) -> str:
    return str(name).replace("_", " ").strip().title()


def _matched_patterns_text(result) -> str:
    matched = tuple(getattr(result, "matched_patterns", ()) or ())
    return ", ".join(_humanize_pattern(name) for name in matched) if matched else "none"


def _format_ranked_result(index: int, result) -> str:
    specificity = cross_category_specificity(result)
    evidence = {
        "high_specificity": "category-enriched evidence",
        "shared": "shared/nonspecific evidence",
        "below": "below category threshold",
        "unavailable": "cross-category calibration unavailable",
    }[specificity]
    return (
        f"{index}. {_clean_label(result.model_label)} | "
        f"score={result.score:.6f} | threshold={result.threshold:.6f} | "
        f"margin={result.margin:.6f} | decision={_clean_decision(result.decision)} | "
        f"cross-category interpretation={evidence} | "
        f"patterns={_matched_patterns_text(result)}"
    )


def format_all_model_results(results) -> str:
    valid_results = [result for result in results if result.valid]
    if not valid_results:
        return "Invalid SMILES. Check syntax and try again."

    product_results = [result for result in valid_results if get_model_role(result.model_id) == PRODUCT_USE_ROLE]
    auxiliary_results = [result for result in valid_results if get_model_role(result.model_id) == AUXILIARY_HAZARD_ROLE]
    rank = {"high_specificity": 2, "shared": 1, "below": 0, "unavailable": -1}
    product_results.sort(
        key=lambda item: (rank[cross_category_specificity(item)], item.score, item.margin),
        reverse=True,
    )
    auxiliary_results.sort(key=lambda item: (item.score, item.margin), reverse=True)

    representative = choose_representative_product_result(product_results)
    positive_products = [result for result in product_results if result.score >= result.threshold]
    if len(positive_products) > 1:
        overlap_text = f"{len(positive_products)} product-use categories are threshold-positive; review overlap/ambiguity."
    else:
        overlap_text = f"{len(positive_products)} product-use category is threshold-positive."

    representative_text = (
        f"Representative product-use evidence: {_clean_label(representative.model_label)}"
        if representative is not None
        else "Representative product-use evidence: unresolved; no single category-enriched signal."
    )
    lines = [
        representative_text,
        "Raw scores and margins are not calibrated probabilities or comparable cross-model distances.",
        overlap_text,
        "",
        "Product-use category scores:",
    ]
    for index, result in enumerate(product_results, start=1):
        lines.append(_format_ranked_result(index, result))

    lines.extend(["", "Auxiliary hazard signal:"])
    if auxiliary_results:
        for index, result in enumerate(auxiliary_results, start=1):
            lines.append(_format_ranked_result(index, result))
    else:
        lines.append("none")
    return "\n".join(lines)


class AlgorithmScoringApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Chemical Category Scorer")
        self.geometry("1080x860")
        self.minsize(820, 560)
        self._image_refs: list[tk.PhotoImage] = []
        self._logo_photo: tk.PhotoImage | None = None
        self._build_widgets()
        self._refresh_model_choices()

    def _build_widgets(self) -> None:
        outer = ttk.Frame(self)
        outer.pack(fill="both", expand=True)

        self.main_canvas = tk.Canvas(outer, highlightthickness=0)
        self.main_canvas.pack(side="left", fill="both", expand=True)
        self.main_scrollbar = ttk.Scrollbar(outer, orient="vertical", command=self.main_canvas.yview)
        self.main_scrollbar.pack(side="right", fill="y")
        self.main_canvas.configure(yscrollcommand=self.main_scrollbar.set)
        self._scroll_enabled = False

        root = ttk.Frame(self.main_canvas, padding=12)
        self._canvas_window = self.main_canvas.create_window((0, 0), window=root, anchor="nw")
        root.bind("<Configure>", self._update_scrollregion)
        self.main_canvas.bind("<Configure>", self._resize_scroll_frame)
        self.main_canvas.bind_all("<MouseWheel>", self._on_mousewheel)

        header = ttk.Frame(root)
        header.pack(fill="x", pady=(0, 12))
        self._load_logo(header)
        header_text_block = ttk.Frame(header)
        header_text_block.pack(side="left", fill="x", expand=True, padx=(12, 0))
        ttk.Label(header_text_block, text="Chemical category scorers", font=("Segoe UI", 15, "bold")).pack(anchor="w")
        ttk.Label(
            header_text_block,
            text="GitHub: https://github.com/phdgil/chemical-category-scorer",
            foreground="#555555",
        ).pack(anchor="w")

        selector = ttk.LabelFrame(root, text="Model selection", padding=12)
        selector.pack(fill="x", pady=(0, 12))
        selector.columnconfigure(1, weight=1)

        ttk.Label(selector, text="Scoring mode").grid(row=0, column=0, sticky="w")
        self.model_choice = tk.StringVar()
        self.selected_model_id = DEFAULT_MODEL_ID
        self.model_box = ttk.Combobox(selector, textvariable=self.model_choice, state="readonly")
        self.model_box.grid(row=0, column=1, sticky="ew", padx=(8, 8))
        self.model_box.bind("<<ComboboxSelected>>", self._sync_selected_model)
        ttk.Button(selector, text="Refresh models", command=self._refresh_model_choices).grid(row=0, column=2)

        single = ttk.LabelFrame(root, text="Single SMILES scoring", padding=12)
        single.pack(fill="x", pady=(0, 12))
        single.columnconfigure(1, weight=1)

        ttk.Label(single, text="SMILES").grid(row=0, column=0, sticky="w")
        self.single_smiles = tk.Text(single, height=4, width=100)
        self.single_smiles.grid(row=1, column=0, columnspan=4, sticky="ew", pady=(4, 8))

        ttk.Button(single, text="Score molecule", command=self.score_single).grid(row=2, column=0, sticky="w")
        ttk.Button(single, text="Clear", command=self.clear_single).grid(row=2, column=1, sticky="w", padx=(8, 0))

        self.single_result = tk.Text(
            single,
            height=18,
            wrap="word",
            bg="#ffffff",
            fg="#222222",
            relief="solid",
            borderwidth=1,
            padx=8,
            pady=8,
        )
        self.single_result.grid(row=3, column=0, columnspan=4, sticky="ew", pady=(10, 0))
        self.single_result.configure(state="disabled", cursor="arrow")
        self._set_single_result("Enter a SMILES string and click 'Score molecule'.")

        visuals = ttk.LabelFrame(root, text="Molecule and matched patterns", padding=12)
        visuals.pack(fill="both", pady=(0, 12))
        visuals.columnconfigure(1, weight=1)

        self.molecule_image_label = ttk.Label(visuals, text="Molecule image will appear after scoring.")
        self.molecule_image_label.grid(row=0, column=0, rowspan=2, sticky="nw")

        self.pattern_summary = tk.Text(
            visuals,
            height=4,
            wrap="word",
            bg="#ffffff",
            fg="#222222",
            relief="solid",
            borderwidth=1,
            padx=8,
            pady=8,
        )
        self.pattern_summary.grid(row=0, column=1, sticky="ew", padx=(12, 0))
        self.pattern_summary.configure(state="disabled", cursor="arrow")
        self._set_pattern_summary("Matched structural patterns will appear here when the selected scorer uses them.")
        self.pattern_gallery = ttk.Frame(visuals)
        self.pattern_gallery.grid(row=1, column=1, sticky="nw", padx=(12, 0), pady=(8, 0))

        batch = ttk.LabelFrame(root, text="Batch CSV scoring", padding=12)
        batch.pack(fill="x", pady=(0, 12))
        batch.columnconfigure(1, weight=1)

        self.input_path = tk.StringVar()
        self.output_path = tk.StringVar(value="")
        self.smiles_column = tk.StringVar(value="SMILES")

        ttk.Label(batch, text="Input CSV").grid(row=0, column=0, sticky="w")
        ttk.Entry(batch, textvariable=self.input_path).grid(row=0, column=1, sticky="ew", padx=(8, 8))
        ttk.Button(batch, text="Browse", command=self.browse_input).grid(row=0, column=2)

        ttk.Label(batch, text="Output CSV").grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Entry(batch, textvariable=self.output_path).grid(row=1, column=1, sticky="ew", padx=(8, 8), pady=(8, 0))
        ttk.Button(batch, text="Save as", command=self.browse_output).grid(row=1, column=2, pady=(8, 0))

        ttk.Label(batch, text="SMILES column").grid(row=2, column=0, sticky="w", pady=(8, 0))
        ttk.Entry(batch, textvariable=self.smiles_column).grid(row=2, column=1, sticky="w", padx=(8, 8), pady=(8, 0))
        ttk.Button(batch, text="Score batch", command=self.score_batch).grid(row=2, column=2, pady=(8, 0))

        notes = ttk.LabelFrame(root, text="Notes", padding=12)
        notes.pack(fill="both", expand=True)
        notes_text = (
            "- Core workflow: paste one SMILES for screening or score a CSV file in batch mode.\n"
            "- The model menu shows the public release models; All models is the default single-molecule view.\n"
            "- All models separates category-enriched evidence from shared evidence and keeps endocrine disruption as an auxiliary hazard signal.\n"
            "- A representative product-use result is reported only when exactly one score reaches its high-specificity cross-category threshold.\n"
            "- Raw scores and margins are screening values, not calibrated probabilities or comparable cross-model distances.\n"
            "- Output CSV starts blank on purpose; after you choose an input file, the app suggests an output name automatically.\n"
            "- The molecule image and matched structural patterns are shown whenever pattern matches exist."
        )
        ttk.Label(notes, text=notes_text, justify="left").pack(anchor="w")

    def _update_scrollregion(self, _event=None) -> None:
        bbox = self.main_canvas.bbox("all")
        if not bbox:
            self._scroll_enabled = False
            return
        self.main_canvas.coords(self._canvas_window, 0, 0)
        self.main_canvas.configure(scrollregion=bbox)
        content_height = max(0, bbox[3] - bbox[1])
        viewport_height = max(1, self.main_canvas.winfo_height())
        self._scroll_enabled = content_height > (viewport_height + 4)
        if self._scroll_enabled:
            if not self.main_scrollbar.winfo_ismapped():
                self.main_scrollbar.pack(side="right", fill="y")
        else:
            self.main_canvas.yview_moveto(0.0)
            if self.main_scrollbar.winfo_ismapped():
                self.main_scrollbar.pack_forget()

    def _resize_scroll_frame(self, event) -> None:
        self.main_canvas.itemconfigure(self._canvas_window, width=event.width)
        self._update_scrollregion()

    def _on_mousewheel(self, event) -> str | None:
        if not self.main_canvas.winfo_exists() or not self._scroll_enabled:
            return "break"
        delta = int(-1 * (event.delta / 120))
        if delta == 0:
            return "break"
        first, last = self.main_canvas.yview()
        if delta < 0 and first <= 0.0:
            return "break"
        if delta > 0 and last >= 1.0:
            return "break"
        self.main_canvas.yview_scroll(delta, "units")
        return "break"

    def _load_logo(self, parent: ttk.Frame) -> None:
        if not LOGO_PATH.exists():
            return
        try:
            original = tk.PhotoImage(file=str(LOGO_PATH), master=self)
            self.iconphoto(True, original)
            max_dim = 96
            factor = max(1, (max(original.width(), original.height()) + max_dim - 1) // max_dim)
            self._logo_photo = original.subsample(factor, factor) if factor > 1 else original
            ttk.Label(parent, image=self._logo_photo).pack(side="left", anchor="nw")
        except Exception:
            self._logo_photo = None

    def _refresh_model_choices(self) -> None:
        refresh_model_registry()
        previous_model_id = self._get_selected_model_id() if self.model_choice.get().strip() else ALL_MODELS_SENTINEL
        self.model_display_to_id = {"All models (screening overview)": ALL_MODELS_SENTINEL}
        for model in list_models(public_only=True):
            display = _clean_label(model["label"])
            if model.get("role") == AUXILIARY_HAZARD_ROLE:
                display = f"{display} (auxiliary hazard)"
            self.model_display_to_id[display] = model["model_id"]
        values = list(self.model_display_to_id.keys())
        self.model_box["values"] = values
        selected_model_id = previous_model_id if previous_model_id in self.model_display_to_id.values() else ALL_MODELS_SENTINEL
        default_display = next((label for label, model_id in self.model_display_to_id.items() if model_id == selected_model_id), values[0])
        self.model_choice.set(default_display)
        self.model_box.set(default_display)
        self.selected_model_id = self.model_display_to_id[default_display]

    def _sync_selected_model(self, _event=None) -> None:
        display = self.model_box.get().strip() or self.model_choice.get().strip()
        mapping = getattr(self, "model_display_to_id", {})
        if display in mapping:
            self.model_choice.set(display)
            self.selected_model_id = mapping[display]

    def _get_selected_model_id(self) -> str:
        self._sync_selected_model()
        display = self.model_box.get().strip() or self.model_choice.get().strip()
        if display in getattr(self, "model_display_to_id", {}):
            return self.model_display_to_id[display]
        return self.selected_model_id if self.selected_model_id in getattr(self, "model_display_to_id", {}).values() else DEFAULT_MODEL_ID

    def _set_readonly_text(self, widget: tk.Text, text: str) -> None:
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", text)
        widget.configure(state="disabled")

    def _set_single_result(self, text: str) -> None:
        self._set_readonly_text(self.single_result, text)

    def _set_pattern_summary(self, text: str) -> None:
        self._set_readonly_text(self.pattern_summary, text)

    def _photo_from_png(self, png_bytes: bytes | None) -> tk.PhotoImage | None:
        if not png_bytes:
            return None
        encoded = base64.b64encode(png_bytes).decode("ascii")
        photo = tk.PhotoImage(data=encoded, master=self)
        self._image_refs.append(photo)
        return photo

    def _clear_visuals(self) -> None:
        self._image_refs = []
        self.molecule_image_label.configure(image="", text="Molecule image will appear after scoring.")
        self._set_pattern_summary("Matched structural patterns will appear here when the selected scorer uses them.")
        for child in self.pattern_gallery.winfo_children():
            child.destroy()

    def _render_visuals(self, smiles: str, result) -> None:
        self._clear_visuals()
        mol_png = render_molecule_png(smiles)
        mol_photo = self._photo_from_png(mol_png)
        if mol_photo is not None:
            self.molecule_image_label.configure(image=mol_photo, text="")
        matched = list(getattr(result, "matched_patterns", ()) or ())
        if not matched:
            self._set_pattern_summary(f"No structural pattern match was found for {_clean_label(result.model_label)}.")
            return
        self._set_pattern_summary(
            f"Matched structural patterns for {_clean_label(result.model_label)}: "
            + ", ".join(_humanize_pattern(name) for name in matched)
        )
        for index, pattern_name in enumerate(matched[:6]):
            frame = ttk.Frame(self.pattern_gallery, padding=(0, 0, 12, 12))
            frame.grid(row=index // 3, column=index % 3, sticky="nw")
            ttk.Label(frame, text=_humanize_pattern(pattern_name)).pack(anchor="w")
            png = render_pattern_match_png(smiles, result.model_id, pattern_name)
            photo = self._photo_from_png(png)
            if photo is not None:
                ttk.Label(frame, image=photo).pack(anchor="w", pady=(4, 0))
            else:
                ttk.Label(frame, text="Pattern image unavailable").pack(anchor="w", pady=(4, 0))

    def browse_input(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if path:
            self.input_path.set(path)
            if not self.output_path.get().strip():
                input_path = Path(path)
                self.output_path.set(str(input_path.with_name(f"{input_path.stem}_scored.csv")))

    def browse_output(self) -> None:
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if path:
            self.output_path.set(path)

    def clear_single(self) -> None:
        self.single_smiles.delete("1.0", "end")
        self._set_single_result("Enter a SMILES string and click 'Score molecule'.")
        self._clear_visuals()

    def score_single(self) -> None:
        smiles = self.single_smiles.get("1.0", "end").strip()
        if not smiles:
            self._set_single_result("Enter a SMILES string first.")
            self._clear_visuals()
            return
        selected = self._get_selected_model_id()
        if selected == ALL_MODELS_SENTINEL:
            results = score_smiles_all(smiles)
            best = choose_best_product_result(results) or choose_best_result(results)
            if not best:
                self._set_single_result("Invalid SMILES. Check syntax and try again.")
                self._clear_visuals()
                return
            self._set_single_result(format_all_model_results(results))
            self._render_visuals(smiles, best)
            return

        result = score_smiles(smiles, selected)
        if not result.valid:
            self._set_single_result("Invalid SMILES. Check syntax and try again.")
            self._clear_visuals()
            return
        matched_text = ", ".join(_humanize_pattern(name) for name in result.matched_patterns) if result.matched_patterns else "none"
        self._set_single_result(
            "\n".join(
                [
                    f"Model: {_clean_label(result.model_label)}",
                    f"Decision: {_clean_decision(result.decision)}",
                    f"Final score: {result.score:.6f}",
                    f"Threshold: {result.threshold:.6f}",
                    f"Margin: {result.margin:.6f}",
                    f"Matched structural patterns: {matched_text}",
                ]
            )
        )
        self._render_visuals(smiles, result)

    def score_batch(self) -> None:
        input_csv = self.input_path.get().strip()
        output_csv = self.output_path.get().strip()
        smiles_column = self.smiles_column.get().strip() or None
        if not input_csv or not output_csv:
            messagebox.showerror("Missing path", "Set both input and output CSV paths.")
            return
        selected = self._get_selected_model_id()
        model_ids = [selected] if selected != ALL_MODELS_SENTINEL else [ALL_MODELS_SENTINEL]
        try:
            out_path = score_csv(input_csv, output_csv, smiles_column, model_ids)
        except Exception as exc:
            messagebox.showerror("Batch scoring failed", str(exc))
            return
        messagebox.showinfo("Batch scoring complete", f"Saved scored CSV to:\n{out_path}")


def run_self_test() -> None:
    refresh_model_registry()
    models = list_models(public_only=True)
    print("SELF_TEST_MODELS=" + ",".join(model["model_id"] for model in models))
    assert len(models) == 4
    assert "final_endocrine_disruptors" not in {model["model_id"] for model in models}
    assert sum(1 for model in models if model["role"] == PRODUCT_USE_ROLE) == 3
    assert [model["model_id"] for model in models if model["role"] == AUXILIARY_HAZARD_ROLE] == ["han_endocrine_disruptors"]

    single = score_smiles("CCO", DEFAULT_MODEL_ID)
    print(f"SELF_TEST_DEFAULT={single.model_id}:{single.score:.6f}:{single.decision}")

    ranked = score_smiles_all("CCO")
    assert len(ranked) == 4
    assert "final_endocrine_disruptors" not in {result.model_id for result in ranked}
    best = choose_best_result(ranked)
    if best:
        print(f"SELF_TEST_BEST={best.model_id}:{best.score:.6f}:{best.margin:.6f}:{len(best.matched_patterns)}")
    formatted_ranked = format_all_model_results(ranked)
    assert "Representative product-use evidence:" in formatted_ranked
    assert "Raw scores and margins are not calibrated probabilities" in formatted_ranked
    assert "Product-use category scores:" in formatted_ranked
    assert "Auxiliary hazard signal:" in formatted_ranked
    assert "threshold=" in formatted_ranked
    assert "patterns=" in formatted_ranked

    assert render_molecule_png("CCO") is not None
    print("SELF_TEST_IMAGE=True")

    app = AlgorithmScoringApp()
    app.withdraw()
    app.update_idletasks()
    default_single = app.single_result.get("1.0", "end").strip()
    default_pattern = app.pattern_summary.get("1.0", "end").strip()
    assert default_single == "Enter a SMILES string and click 'Score molecule'."
    assert default_pattern == "Matched structural patterns will appear here when the selected scorer uses them."
    assert app._get_selected_model_id() == ALL_MODELS_SENTINEL
    app.model_box.set("Flavor and fragrance category score")
    app._sync_selected_model()
    app.single_smiles.insert("1.0", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC")
    app.score_single()
    app.update_idletasks()
    scored_single = app.single_result.get("1.0", "end").strip()
    scored_pattern = app.pattern_summary.get("1.0", "end").strip()
    assert "Model: Flavor and fragrance category score" in scored_single
    assert "Decision:" in scored_single
    assert "Threshold:" in scored_single
    assert "Matched structural patterns:" in scored_single
    assert scored_pattern.startswith("Matched structural patterns for Flavor and fragrance category score") or scored_pattern.startswith("No structural pattern match was found for Flavor and fragrance category score")
    print("SELF_TEST_UI_TEXT=True")
    app.destroy()

    with tempfile.TemporaryDirectory() as tmp:
        input_path = Path(tmp) / "input.csv"
        output_path = Path(tmp) / "output.csv"
        with input_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["SMILES", "name"])
            writer.writeheader()
            writer.writerow({"SMILES": "CCO", "name": "ethanol"})
            writer.writerow({"SMILES": "CC(=O)Oc1ccccc1C(=O)O", "name": "aspirin"})
        score_csv(input_path, output_path, "SMILES", [ALL_MODELS_SENTINEL])
        print(f"SELF_TEST_BATCH={output_path.exists()}")
        print(output_path.read_text(encoding="utf-8").strip())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Desktop app for chemical category scoring.")
    parser.add_argument("--self-test", action="store_true", help="Run a headless smoke test.")
    parser.add_argument("--list-models", action="store_true", help="Print available model ids and labels.")
    parser.add_argument("--score", help="Score a single SMILES string with one model.")
    parser.add_argument("--model-id", default=DEFAULT_MODEL_ID, help="Model id for --score or --batch-*.")
    parser.add_argument("--score-all", help="Score a single SMILES string across all models.")
    parser.add_argument("--batch-in", help="Input CSV for batch scoring.")
    parser.add_argument("--batch-out", help="Output CSV for batch scoring.")
    parser.add_argument("--smiles-column", default="SMILES", help="SMILES column name for batch scoring.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.list_models:
        refresh_model_registry()
        for model in list_models(public_only=True):
            print(f"{model['model_id']}\t{_clean_label(model['label'])}")
        return
    if args.score_all:
        for result in score_smiles_all(args.score_all):
            print(result)
        return
    if args.score:
        print(score_smiles(args.score, args.model_id))
        return
    if args.batch_in and args.batch_out:
        model_ids = [args.model_id] if args.model_id != ALL_MODELS_SENTINEL else [ALL_MODELS_SENTINEL]
        print(score_csv(args.batch_in, args.batch_out, args.smiles_column, model_ids))
        return
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    app = AlgorithmScoringApp()
    app.mainloop()


if __name__ == "__main__":
    main()
