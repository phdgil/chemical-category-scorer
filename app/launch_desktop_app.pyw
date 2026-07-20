from __future__ import annotations

import os
import sys
import threading
import traceback
from pathlib import Path
import tkinter as tk
from tkinter import messagebox, ttk

APP_DIR = Path(__file__).resolve().parent
LOG_PATH = Path(os.environ.get("TEMP", str(APP_DIR))) / "chemical_category_scorer_launcher.log"


def _show_error(message: str) -> None:
    try:
        root = tk.Tk()
        root.withdraw()
        messagebox.showerror("Chemical Category Scorer", message)
        root.destroy()
    except Exception:
        pass


def _make_splash() -> tuple[tk.Tk, tk.StringVar]:
    splash = tk.Tk()
    splash.title("Chemical Category Scorer")
    splash.geometry("440x130")
    splash.resizable(False, False)
    frame = ttk.Frame(splash, padding=16)
    frame.pack(fill="both", expand=True)
    ttk.Label(frame, text="Chemical Category Scorer", font=("Segoe UI", 12, "bold")).pack(anchor="w")
    status = tk.StringVar(value="Loading models and starting the desktop app...")
    ttk.Label(frame, textvariable=status, padding=(0, 10, 0, 0)).pack(anchor="w")
    ttk.Label(frame, text="Please wait. The window stays responsive while the app loads.", foreground="#555555").pack(anchor="w")
    splash.update_idletasks()
    try:
        width = splash.winfo_width()
        height = splash.winfo_height()
        x = (splash.winfo_screenwidth() - width) // 2
        y = (splash.winfo_screenheight() - height) // 2
        splash.geometry(f"{width}x{height}+{x}+{y}")
    except Exception:
        pass
    return splash, status


def _background_import(state: dict[str, object]) -> None:
    try:
        import desktop_app  # type: ignore
        state["module"] = desktop_app
        state["done"] = True
    except Exception:
        state["error"] = traceback.format_exc()
        state["done"] = True


def main() -> None:
    splash = None
    try:
        os.chdir(APP_DIR)
        sys.path.insert(0, str(APP_DIR))
        splash, status = _make_splash()
        state: dict[str, object] = {"done": False, "module": None, "error": None}
        worker = threading.Thread(target=_background_import, args=(state,), daemon=True)
        worker.start()

        ticks = {"count": 0}

        def poll() -> None:
            ticks["count"] += 1
            dots = "." * (ticks["count"] % 4)
            status.set(f"Loading models and starting the desktop app{dots}")
            if state["done"]:
                splash.quit()
                return
            splash.after(150, poll)

        splash.after(150, poll)
        splash.mainloop()

        if state.get("error"):
            raise RuntimeError(str(state["error"]))

        module = state.get("module")
        if module is None:
            raise RuntimeError("Desktop app import finished without a module result.")

        status.set("Opening the main window...")
        splash.update_idletasks()
        app = module.AlgorithmScoringApp()
        app.update_idletasks()
        splash.destroy()
        splash = None
        app.mainloop()
    except Exception:
        if splash is not None:
            try:
                splash.destroy()
            except Exception:
                pass
        LOG_PATH.write_text(traceback.format_exc(), encoding="utf-8")
        _show_error(f"The desktop app could not start.\n\nSee log:\n{LOG_PATH}")
        raise


if __name__ == "__main__":
    main()
