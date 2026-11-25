# -*- coding: utf-8 -*-
"""Ce 3d XPS Fitter v2: a responsive, cancellation-aware Tk interface."""
from __future__ import annotations

import argparse
import json
import math
import os
import queue
import re
import sys
import threading
import traceback
from datetime import datetime
from pathlib import Path
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, ttk

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
except ImportError:
    DND_FILES = None
    TkinterDnD = None


VERSION = "2.0.1"
SUPPORTED_SUFFIXES = {".xls", ".xlsx", ".vms", ".vamas", ".csv", ".txt", ".dat", ".xy"}
FILE_TYPES = [
    ("XPS data (automatic format detection)", "*.xls *.xlsx *.vms *.vamas *.csv *.txt *.dat *.xy"),
    ("Excel workbooks", "*.xls *.xlsx"),
    ("VAMAS spectra", "*.vms *.vamas"),
    ("Text data", "*.csv *.txt *.dat *.xy"),
    ("All files", "*.*"),
]


def _saved_result_records(result_dir: Path) -> list[tuple[Path, dict]]:
    """Prefer the organized result once per run, with legacy-layout fallback."""
    candidates = [result_dir / "fit_result.json", result_dir / "details" / "fit_result.json"]
    candidates.extend(sorted(result_dir.glob("*/fit_result.json")))
    candidates.extend(sorted(result_dir.glob("*/*/fit_result.json")))
    canonical: dict[Path, Path] = {}
    for path in candidates:
        if not path.is_file():
            continue
        run = path.parent.parent if path.parent.name == "details" else path.parent
        organized = run / "details" / "fit_result.json"
        canonical[run] = organized if organized.is_file() else path
    records = []
    for path in canonical.values():
        try:
            item = json.loads(path.read_text(encoding="utf-8-sig"))
            if isinstance(item, dict):
                records.append((path, item))
        except (OSError, ValueError):
            continue
    return records


def inspect_result(result_dir: Path, requested_count: int) -> tuple[str, str]:
    """Read saved statuses; completion of a worker is not proof of convergence."""
    summary_path = result_dir / "details" / "batch_summary.json"
    if not summary_path.is_file():
        summary_path = result_dir / "batch_summary.json"
    summary = None
    if summary_path.is_file():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8-sig"))
        except (OSError, ValueError):
            pass
    statuses: list[str] = []
    if isinstance(summary, dict):
        entries = summary.get("results", [])
        if isinstance(entries, list):
            statuses = [str(item.get("status", "unknown")) for item in entries if isinstance(item, dict)]
        if summary.get("status") == "cancelled":
            return "cancelled", "Stopped; available results saved. Convergence was not established."
        if summary.get("status") == "completed_with_errors":
            return "failed", "Batch finished with failed files; review the run log and batch summary."
    if not statuses:
        statuses = [str(item.get("status", "unknown")) for _, item in _saved_result_records(result_dir)]
    if "cancelled" in statuses:
        return "cancelled", "Stopped; available results saved. Convergence was not established."
    if "failed" in statuses:
        return "failed", "Run finished with failed files; review the results and log."
    if "budget_exhausted" in statuses:
        return "incomplete", "Run ended before the engineering convergence criteria were met; review the saved status."
    complete = len(statuses) == requested_count
    if complete and all(status == "converged" for status in statuses):
        return "converged", "Engineering convergence criteria met. Results saved."
    if complete and all(status in {"quick_complete", "converged"} for status in statuses):
        return "quick_complete", "Quick fit complete. Results saved; review diagnostics to assess reliability."
    return "unknown", "Run returned; review the actual status in the result files."


class Ce3dXpsFitterApp:
    def __init__(self, root: tk.Tk | None = None) -> None:
        self.dnd_available = False
        if root is not None:
            self.root = root
        elif TkinterDnD is not None:
            try:
                self.root = TkinterDnD.Tk()
                self.dnd_available = True
            except Exception:
                self.root = tk.Tk()
        else:
            self.root = tk.Tk()
        self.root.title(f"Ce 3d XPS Fitter v{VERSION} - English")
        self.root.geometry("900x800")
        self.root.minsize(820, 680)
        self.root.protocol("WM_DELETE_WINDOW", self._request_close)

        self.model = tk.StringVar(self.root, value="LA")
        self.satellites = tk.BooleanVar(self.root, value=False)
        self.mode = tk.StringVar(self.root, value="quick")
        self.output_dir = tk.StringVar(self.root, value=str(Path.home() / "Documents" / "Ce3d_XPS_Fit_Outputs"))
        self.walkers = tk.StringVar(self.root, value="80")
        self.quick_steps = tk.StringVar(self.root, value="8000")
        self.burn = tk.StringVar(self.root, value="4000")
        self.seed = tk.StringVar(self.root, value="42")
        self.vamas_block = tk.StringVar(self.root, value="")
        self.vamas_ordinate = tk.StringVar(self.root, value="")
        self.kinetic_offset = tk.StringVar(self.root, value="")
        self.satellite_gamma = tk.StringVar(self.root, value="0.8")
        self.background_mode = tk.StringVar(self.root, value="fixed_shirley")
        self.status = tk.StringVar(self.root, value="Ready")
        self.step_status = tk.StringVar(self.root, value="Defaults: 80 walkers · 8000 steps · 4000 burn-in steps")
        self.file_count = tk.StringVar(self.root, value="No files selected")
        self.input_files: list[Path] = []
        self.events: queue.Queue[tuple[str, object]] = queue.Queue()
        self.worker: threading.Thread | None = None
        self.stop_event = threading.Event()
        self.active = False
        self.close_when_finished = False
        self.advanced_visible = False
        self.last_output: Path | None = None
        self.session_log: Path | None = None
        self._log_lock = threading.Lock()
        self._controls: list[tk.Widget] = []
        self._inspecting_vamas = False
        self._build_ui()
        self.root.after(100, self._drain_events)

    def mainloop(self) -> None:
        self.root.mainloop()

    def _control(self, widget: tk.Widget) -> tk.Widget:
        self._controls.append(widget)
        return widget

    def _build_ui(self) -> None:
        self.root.option_add("*Font", "{Segoe UI} 9")
        ttk.Style(self.root).configure(".", font=("Segoe UI", 9))
        frame = ttk.Frame(self.root, padding=(16, 12))
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(8, weight=1)
        ttk.Label(frame, text=f"Ce 3d XPS Fitter  {VERSION}", font=("Segoe UI", 15, "bold")).grid(
            row=0, column=0, sticky="w", pady=(0, 8))

        choices = ttk.Frame(frame, padding=(10, 6))
        self.selection_frame = choices
        choices.grid(row=1, column=0, sticky="ew")
        choices.columnconfigure(1, weight=1, uniform="choices")
        choices.columnconfigure(2, weight=1, uniform="choices")
        for row, label, variable, options in (
            (0, "Peak model", self.model, (("LA (default)", "LA"), ("Voigt", "Voigt"))),
            (1, "Satellite peaks", self.satellites, (("Exclude (default)", False), ("Include", True))),
            (2, "Sampling mode", self.mode, (("Quick fit (8000 steps)", "quick"), ("Convergence validation", "convergence"))),
        ):
            ttk.Label(choices, text=label, width=15).grid(row=row, column=0, sticky="w", pady=5)
            for column, (text, value) in enumerate(options, start=1):
                self._control(ttk.Radiobutton(choices, text=text, variable=variable, value=value,
                                              command=self._sync_mode)).grid(row=row, column=column, sticky="w", padx=8)
        ttk.Label(frame, text="Gaussian-convolved LA-type: fixed alpha=0.9, beta=2.",
                  foreground="#555555", wraplength=840).grid(row=2, column=0, sticky="w", pady=(4, 5))

        file_box = ttk.LabelFrame(frame, text="Input files · automatic format detection; files fitted sequentially", padding=8)
        self.input_frame = file_box
        file_box.grid(row=3, column=0, sticky="ew", pady=4)
        file_box.columnconfigure(0, weight=1)
        file_box.rowconfigure(0, weight=1)
        self.file_list = tk.Listbox(file_box, height=4, selectmode=tk.EXTENDED, activestyle="none", exportselection=False)
        self.file_list.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(file_box, orient="vertical", command=self.file_list.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        horizontal = ttk.Scrollbar(file_box, orient="horizontal", command=self.file_list.xview)
        horizontal.grid(row=1, column=0, sticky="ew")
        self.file_list.configure(yscrollcommand=scroll.set, xscrollcommand=horizontal.set)
        buttons = ttk.Frame(file_box)
        buttons.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(6, 0))
        for label, command in (("Add files...", self._browse_files), ("Remove selected", self._remove_files), ("Clear", self._clear_files)):
            self._control(ttk.Button(buttons, text=label, command=command)).pack(side="left", padx=(0, 8))
        self.vamas_select_button = self._control(ttk.Button(buttons, text="Select VAMAS block...", command=self._choose_vamas))
        self.vamas_select_button.pack(side="left", padx=(0, 8))
        ttk.Label(buttons, textvariable=self.file_count).pack(side="right")
        if self.dnd_available:
            try:
                self.file_list.drop_target_register(DND_FILES)
                self.file_list.dnd_bind("<<Drop>>", self._on_drop)
            except Exception:
                self.dnd_available = False

        out_box = ttk.LabelFrame(frame, text="Output folder", padding=8)
        self.output_frame = out_box
        out_box.grid(row=4, column=0, sticky="ew", pady=4)
        out_box.columnconfigure(0, weight=1)
        self._control(ttk.Entry(out_box, textvariable=self.output_dir)).grid(row=0, column=0, sticky="ew", padx=(0, 8))
        self._control(ttk.Button(out_box, text="Browse...", command=self._browse_output)).grid(row=0, column=1)

        advanced_holder = ttk.Frame(frame)
        advanced_holder.grid(row=5, column=0, sticky="ew", pady=3)
        advanced_holder.columnconfigure(0, weight=1)
        self._control(ttk.Button(advanced_holder, text="▸ Advanced settings", command=self._toggle_advanced))
        self.advanced_button = self._controls[-1]
        self.advanced_button.grid(row=0, column=0, sticky="w")
        self.advanced_frame = ttk.LabelFrame(advanced_holder, text="Sampling and format options", padding=8)
        self.advanced_frame.grid(row=1, column=0, sticky="ew", pady=(4, 0))
        for col in (1, 3, 5, 7):
            self.advanced_frame.columnconfigure(col, weight=1)
        fields = (
            ("Walkers", self.walkers), ("Quick steps", self.quick_steps),
            ("Burn-in steps", self.burn), ("Seed", self.seed),
            ("VAMAS block", self.vamas_block), ("Signal index", self.vamas_ordinate),
            ("KE offset / eV", self.kinetic_offset), ("Satellite gamma / eV", self.satellite_gamma),
        )
        for index, (label, variable) in enumerate(fields):
            row, col = divmod(index, 4)
            ttk.Label(self.advanced_frame, text=label).grid(row=row, column=col * 2, sticky="w", padx=(0, 4), pady=3)
            entry = self._control(ttk.Entry(self.advanced_frame, textvariable=variable, width=8))
            entry.grid(row=row, column=col * 2 + 1, sticky="ew", padx=(0, 8), pady=3)
            if variable is self.quick_steps:
                self.quick_steps_entry = entry
        ttk.Label(self.advanced_frame, text="Steps / burn-in: multiples of 2000. VAMAS indices start at 0; leave blank for automatic selection.",
                  foreground="#555555").grid(row=2, column=0, columnspan=8, sticky="w", pady=(3, 2))
        bg_row = ttk.Frame(self.advanced_frame)
        bg_row.grid(row=3, column=0, columnspan=8, sticky="w")
        ttk.Label(bg_row, text="Background:").pack(side="left")
        for text, value in (("Fixed Shirley (default)", "fixed_shirley"), ("Use background from file", "provided")):
            self._control(ttk.Radiobutton(bg_row, text=text, variable=self.background_mode, value=value)).pack(side="left", padx=6)
        self.advanced_frame.grid_remove()

        run_row = ttk.Frame(frame)
        run_row.grid(row=6, column=0, sticky="ew", pady=(6, 2))
        self.run_button = ttk.Button(run_row, text="Start fit", command=self._start_fit)
        self.run_button.pack(side="left")
        self.stop_button = ttk.Button(run_row, text="Stop and save", command=self._request_stop, state="disabled")
        self.stop_button.pack(side="left", padx=8)
        self.open_button = ttk.Button(run_row, text="Open report", command=self._open_results, state="disabled")
        self.open_button.pack(side="left")
        self.status_label = ttk.Label(run_row, textvariable=self.status, foreground="#555555", wraplength=480)
        self.status_label.pack(side="left", padx=12)
        self.status_label.bind("<Configure>", lambda event: self.status_label.configure(wraplength=max(1, event.width)))

        progress_box = ttk.Frame(frame)
        progress_box.grid(row=7, column=0, sticky="ew", pady=(5, 3))
        progress_box.columnconfigure(0, weight=1)
        self.progress_bar = ttk.Progressbar(progress_box, mode="determinate", maximum=8000)
        self.progress_bar.grid(row=0, column=0, sticky="ew")
        self.step_status_label = ttk.Label(progress_box, textvariable=self.step_status, foreground="#555555", wraplength=780)
        self.step_status_label.grid(row=1, column=0, sticky="ew", pady=(3, 0))
        self.step_status_label.bind("<Configure>", lambda event: self.step_status_label.configure(wraplength=max(1, event.width)))

        log_box = ttk.LabelFrame(frame, text="Run log · updates at the start and every 2000 steps", padding=8)
        self.log_frame = log_box
        log_box.grid(row=8, column=0, sticky="nsew", pady=(5, 0))
        log_box.columnconfigure(0, weight=1)
        log_box.rowconfigure(0, weight=1)
        self.log_text = tk.Text(log_box, height=8, wrap="word", state="disabled", font=("Segoe UI", 9))
        self.log_text.grid(row=0, column=0, sticky="nsew")
        log_scroll = ttk.Scrollbar(log_box, orient="vertical", command=self.log_text.yview)
        log_scroll.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=log_scroll.set)
        self._append_log("Supported: .xls / .xlsx / .vms / .vamas / .csv / .txt / .dat / .xy.")
        self._append_log("Each input file and each run use a separate results folder.")
        if self.dnd_available:
            self._append_log("You can drag data files into the input list above.")

    def _toggle_advanced(self) -> None:
        self.advanced_visible = not self.advanced_visible
        if self.advanced_visible:
            # Keep the input/output/log areas visible when advanced controls open.
            self.root.minsize(820, 800)
            self.advanced_frame.grid()
            self.advanced_button.configure(text="▾ Hide advanced settings")
        else:
            self.root.minsize(820, 680)
            self.advanced_frame.grid_remove()
            self.advanced_button.configure(text="▸ Advanced settings")

    def _sync_mode(self) -> None:
        self.quick_steps_entry.configure(state="normal" if self.mode.get() == "quick" and not self.active else "disabled")
        if not self.active:
            if self.mode.get() == "convergence":
                self.step_status.set("Convergence validation: runs until the criteria are met, with no preset step limit. Stop and save is available.")
            else:
                self.step_status.set("Quick fit: uses the steps in Advanced settings. Completion does not establish convergence.")

    def _browse_files(self) -> None:
        self._add_files(filedialog.askopenfilenames(parent=self.root, title="Select one or more Ce 3d data files", filetypes=FILE_TYPES))

    def _on_drop(self, event) -> None:
        if not self.active:
            self._add_files(self.root.tk.splitlist(event.data))

    def _add_files(self, files) -> None:
        if self.active:
            return
        rejected = []
        for item in files:
            path = Path(str(item)).expanduser().resolve()
            if not path.is_file() or path.suffix.lower() not in SUPPORTED_SUFFIXES:
                rejected.append(str(path))
                continue
            if path not in self.input_files:
                self.input_files.append(path)
                self.file_list.insert(tk.END, str(path))
        self.file_count.set(f"{len(self.input_files)} file(s) selected" if self.input_files else "No files selected")
        if rejected:
            self._append_log("Ignored unsupported formats or missing files:\n" + "\n".join(rejected))

    def _remove_files(self) -> None:
        if self.active:
            return
        for index in reversed(self.file_list.curselection()):
            del self.input_files[index]
            self.file_list.delete(index)
        self.file_count.set(f"{len(self.input_files)} file(s) selected" if self.input_files else "No files selected")

    def _clear_files(self) -> None:
        if not self.active:
            self.input_files.clear()
            self.file_list.delete(0, tk.END)
            self.file_count.set("No files selected")

    def _browse_output(self) -> None:
        dirname = filedialog.askdirectory(parent=self.root, title="Select output folder", initialdir=self.output_dir.get() or str(Path.home()))
        if dirname:
            self.output_dir.set(dirname)

    def _choose_vamas(self) -> None:
        if self.active or self._inspecting_vamas:
            return
        selected = [self.input_files[index] for index in self.file_list.curselection()
                    if self.input_files[index].suffix.lower() in {".vms", ".vamas"}]
        available = [path for path in self.input_files if path.suffix.lower() in {".vms", ".vamas"}]
        if len(selected) == 1:
            path = selected[0]
        elif len(available) == 1:
            path = available[0]
        else:
            messagebox.showinfo("Select a VAMAS file", "Add a VAMAS file, then select the file to inspect in the input list.", parent=self.root)
            return
        self._inspecting_vamas = True
        self.vamas_select_button.configure(state="disabled")
        self.run_button.configure(state="disabled")
        self.status.set("Reading VAMAS block information...")

        def inspect() -> None:
            try:
                from ce3d_vamas import inspect_vamas
                self.events.put(("vamas_inventory", (path, inspect_vamas(path))))
            except Exception as exc:
                self.events.put(("vamas_error", (str(exc), traceback.format_exc())))

        # This read-only helper never owns sampler state or writes result files.
        threading.Thread(target=inspect, name="VamasInventory", daemon=True).start()

    def _show_vamas_choices(self, path: Path, inventory: dict) -> None:
        blocks = [item for item in inventory.get("block_inventory", [])
                  if str(item.get("technique", "")).strip().upper() == "XPS"]
        candidates = set(inventory.get("ce3d_candidate_indices", []))
        if not blocks:
            messagebox.showerror("Cannot select a block", "This VAMAS file contains no selectable XPS blocks.", parent=self.root)
            return
        dialog = tk.Toplevel(self.root)
        dialog.title("Select a VAMAS block and signal channel")
        dialog.geometry("850x440")
        dialog.minsize(820, 440)
        dialog.transient(self.root)
        body = ttk.Frame(dialog, padding=14)
        body.pack(fill="both", expand=True)
        body.columnconfigure(0, weight=1)
        ttk.Label(body, text=path.name, font=("Segoe UI", 11, "bold"), wraplength=750).grid(row=0, column=0, sticky="w")
        ttk.Label(body, text="Select the measured Ce 3d block and intensity channel. These settings apply to all VAMAS files in this run.\nRun files separately if their block or channel layouts differ.",
                  wraplength=750, foreground="#555555").grid(row=1, column=0, sticky="w", pady=(5, 12))
        ttk.Label(body, text="XPS block (Ce 3d candidates marked; indices start at 0)").grid(row=2, column=0, sticky="w")
        block_labels = [f"[{item['index']}] {'Ce 3d candidate · ' if item['index'] in candidates else ''}{item.get('label', '')}  |  {item.get('sample_identifier', '')}"
                        for item in blocks]
        block_box = ttk.Combobox(body, state="readonly", values=block_labels, width=90)
        block_box.grid(row=3, column=0, sticky="ew", pady=(3, 6))
        details = tk.StringVar(dialog, value="Select a block. When several candidates exist, none is selected automatically.")
        ttk.Label(body, textvariable=details, wraplength=750, foreground="#555555").grid(row=4, column=0, sticky="w", pady=(0, 8))
        ttk.Label(body, text="Signal channel (measured intensity; transmission / TXFN is not an intensity signal)").grid(row=5, column=0, sticky="w")
        ordinate_box = ttk.Combobox(body, state="readonly", width=90)
        ordinate_box.grid(row=6, column=0, sticky="ew", pady=(3, 8))
        offset_row = ttk.Frame(body)
        offset_row.grid(row=7, column=0, sticky="w", pady=4)
        ttk.Label(offset_row, text="KE offset / eV:").pack(side="left")
        offset = tk.StringVar(dialog, value=self.kinetic_offset.get())
        offset_entry = ttk.Entry(offset_row, textvariable=offset, width=14)
        offset_entry.pack(side="left")
        ttk.Label(offset_row, text="For kinetic-energy axes only; confirm the instrument convention.", foreground="#555555").pack(side="left", padx=8)
        button_row = ttk.Frame(body)
        button_row.grid(row=8, column=0, sticky="e", pady=(12, 0))
        current_ordinates: list[dict] = []

        def block_changed(_event=None) -> None:
            index = block_box.current()
            if index < 0:
                return
            block = blocks[index]
            details.set(f"{block.get('axis', '')} / {block.get('units', '')} · "
                        f"{block.get('start', '')} → {block.get('end', '')} · step {block.get('step', '')} · "
                        f"{block.get('point_count', '')} points")
            current_ordinates[:] = block.get("ordinates", [])
            ordinate_box.configure(values=[f"[{item['index']}] {item.get('label', '')} / {item.get('unit', '')}"
                                           for item in current_ordinates])
            ordinate_box.set("")
            if len(current_ordinates) == 1:
                ordinate_box.current(0)
            elif self.vamas_block.get().strip() == str(block["index"]):
                for pos, item in enumerate(current_ordinates):
                    if self.vamas_ordinate.get().strip() == str(item["index"]):
                        ordinate_box.current(pos)
                        break
            kinetic = "kinetic" in str(block.get("axis", "")).lower()
            offset_entry.configure(state="normal" if kinetic else "disabled")

        def apply() -> None:
            if block_box.current() < 0 or ordinate_box.current() < 0:
                messagebox.showerror("Incomplete selection", "Select an XPS block and a measured-intensity channel.", parent=dialog)
                return
            block = blocks[block_box.current()]
            ordinate = current_ordinates[ordinate_box.current()]
            signal_label = str(ordinate.get("label", "")).lower()
            if "transmission" in signal_label or "txfn" in signal_label:
                messagebox.showerror("Select an intensity channel", "The selected channel is a transmission function, not measured spectral intensity.", parent=dialog)
                return
            kinetic = "kinetic" in str(block.get("axis", "")).lower()
            try:
                if kinetic:
                    self._number(offset.get(), "KE offset")
            except ValueError as exc:
                messagebox.showerror("Incomplete kinetic-energy conversion", str(exc), parent=dialog)
                return
            self.vamas_block.set(str(block["index"]))
            self.vamas_ordinate.set(str(ordinate["index"]))
            self.kinetic_offset.set(offset.get().strip() if kinetic else "")
            self._append_log(f"VAMAS selection: {path.name}; block [{block['index']}] {block.get('label', '')};"
                             f" channel [{ordinate['index']}] {ordinate.get('label', '')}.")
            if not self.advanced_visible:
                self._toggle_advanced()
            dialog.destroy()

        block_box.bind("<<ComboboxSelected>>", block_changed)
        preselected = self.vamas_block.get().strip()
        for pos, block in enumerate(blocks):
            if str(block["index"]) == preselected or (not preselected and len(candidates) == 1 and block["index"] in candidates):
                block_box.current(pos)
                block_changed()
                break
        ttk.Button(button_row, text="Cancel", command=dialog.destroy).pack(side="left", padx=8)
        ttk.Button(button_row, text="Apply selection", command=apply).pack(side="left")
        dialog.grab_set()

    @staticmethod
    def _integer(text: str, name: str, minimum: int) -> int:
        try:
            value = int(text.strip())
        except ValueError:
            raise ValueError(f"{name} must be an integer.") from None
        if value < minimum:
            raise ValueError(f"{name} must be at least {minimum}.")
        return value

    @staticmethod
    def _number(text: str, name: str, positive: bool = False) -> float:
        try:
            value = float(text.strip())
        except ValueError:
            raise ValueError(f"{name} must be a number.") from None
        if not math.isfinite(value) or (positive and value <= 0):
            raise ValueError(f"{name} must be a finite {'positive number' if positive else 'number'}.")
        return value

    def _snapshot_arguments(self) -> tuple[list[Path], dict]:
        """Called only on the Tk thread; the worker never reads Tk variables."""
        if not self.input_files:
            raise ValueError("Select at least one input file.")
        paths = list(self.input_files)
        if any(not path.is_file() for path in paths):
            raise ValueError("Some input files no longer exist; select them again.")
        output = self.output_dir.get().strip()
        if not output:
            raise ValueError("Select an output folder.")
        mode = self.mode.get()
        walkers = self._integer(self.walkers.get(), "Walkers", 22)
        if walkers % 2:
            raise ValueError("Walkers must be an even number.")
        # In convergence mode this value is not a stopping budget.
        steps = self._integer(self.quick_steps.get(), "Quick-fit steps", 1) if mode == "quick" else 8000
        burn = self._integer(self.burn.get(), "Burn-in steps", 2000)
        if burn % 2000:
            raise ValueError("Burn-in steps must be a multiple of 2000.")
        if mode == "quick" and steps % 2000:
            raise ValueError("Quick-fit steps must be a multiple of 2000.")
        if mode == "quick" and burn >= steps:
            raise ValueError("Burn-in steps must be less than the total steps for a quick fit.")
        kwargs = dict(
            outdir=Path(output).expanduser().resolve(), model=self.model.get(),
            satellites=self.satellites.get(), mode=mode, walkers=walkers, steps=steps,
            burn=burn, seed=self._integer(self.seed.get(), "Random seed", 0),
            background_mode=self.background_mode.get(), max_steps=None,
            satellite_gamma=self._number(self.satellite_gamma.get(), "Satellite gamma", positive=True),
            vamas_block=self._integer(self.vamas_block.get(), "VAMAS block index", 0) if self.vamas_block.get().strip() else None,
            vamas_ordinate=self._integer(self.vamas_ordinate.get(), "Signal channel index", 0) if self.vamas_ordinate.get().strip() else None,
            kinetic_offset_ev=self._number(self.kinetic_offset.get(), "KE offset") if self.kinetic_offset.get().strip() else None,
        )
        if kwargs["seed"] >= 2**32:
            raise ValueError("Random seed must be between 0 and 4294967295.")
        kwargs["outdir"].mkdir(parents=True, exist_ok=True)
        return paths, kwargs

    def _set_running(self, running: bool) -> None:
        self.active = running
        for widget in self._controls:
            widget.configure(state="disabled" if running else "normal")
        self.file_list.configure(state="disabled" if running else "normal")
        self.run_button.configure(state="disabled" if running else "normal")
        self.stop_button.configure(state="normal" if running else "disabled")
        self.open_button.configure(state="normal" if not running and self.last_output is not None else "disabled")
        if not running:
            self.quick_steps_entry.configure(state="normal" if self.mode.get() == "quick" else "disabled")

    def _start_fit(self) -> None:
        if self.active or self._inspecting_vamas or (self.worker is not None and self.worker.is_alive()):
            return
        try:
            paths, kwargs = self._snapshot_arguments()
            log_dir = kwargs["outdir"] / "gui_logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            self.session_log = log_dir / f"gui_{datetime.now():%Y%m%d_%H%M%S_%f}.log"
            self.session_log.write_text(f"Ce 3d XPS Fitter v{VERSION}\n", encoding="utf-8")
        except Exception as exc:
            messagebox.showerror("Invalid input settings", str(exc), parent=self.root)
            return
        self.stop_event = threading.Event()
        self.close_when_finished = False
        self.last_output = None
        self._set_running(True)
        self.step_status_label.configure(font=("Segoe UI", 9), foreground="#555555")
        self.status.set("Starting fit...")
        self.status_label.configure(foreground="#555555")
        self.progress_bar.stop()
        self.progress_bar.configure(value=0)
        if kwargs["mode"] == "convergence":
            self.progress_bar.configure(mode="indeterminate")
            self.progress_bar.start(25)
            self.step_status.set("Convergence validation: waiting for sampling progress; no preset step limit.")
        else:
            self.progress_bar.configure(mode="determinate", maximum=kwargs["steps"])
            self.step_status.set(f"Waiting for sampling · target: {kwargs['steps']} steps")
        self._append_log("\n========== New fit run ==========")
        self._append_log(f"Run log: {self.session_log}")
        # A non-daemon worker allows cancellation/export to finish before exit.
        self.worker = threading.Thread(target=self._fit_worker, args=(paths, kwargs, self.session_log, self.stop_event),
                                       name="Ce3dFitWorker", daemon=False)
        self.worker.start()

    def _fit_worker(self, paths: list[Path], kwargs: dict, log_path: Path, stop_event: threading.Event) -> None:
        def log(message) -> None:
            text = str(message)
            self.events.put(("log", text))
            try:
                with self._log_lock:
                    with log_path.open("a", encoding="utf-8") as handle:
                        handle.write(text + "\n")
            except OSError:
                pass

        def progress(info: dict) -> None:
            self.events.put(("progress", dict(info)))

        try:
            from ce3d_xps_fitter_core import run_fit

            log(f"Model={kwargs['model']}; satellites={'include' if kwargs['satellites'] else 'exclude'}; mode={kwargs['mode']}")
            log(f"Files={len(paths)}; output folder={kwargs['outdir']}")
            if kwargs["mode"] == "convergence":
                log("Convergence validation runs until the criteria are met; no maximum step count is set. Stop requests saving available samples.")
            result = Path(run_fit(paths, **kwargs, log=log, progress=progress, stop_event=stop_event))
            state, text = inspect_result(result, len(paths))
            from ce3d_results import result_headline
            result_lines = []
            for result_path, saved in _saved_result_records(result):
                if not saved.get("posterior_available"):
                    continue
                headline = result_headline(saved)
                quantiles = saved.get("Ce3_percent", {})
                detail = "Posterior mean ± 1 SD (percentage points)"
                if all(key in quantiles for key in ("q025", "q975", "median")):
                    detail += (f" · 95% credible interval {quantiles['q025']:.2f}–{quantiles['q975']:.2f}%"
                               f" · median {quantiles['median']:.2f}%")
                source = saved.get("configuration", {}).get("source")
                run_name = result_path.parent.parent.name if result_path.parent.name == "details" else result_path.parent.name
                sample_name = Path(source).name if source else run_name
                log(f"{sample_name}: {headline}" if len(paths) > 1 else headline)
                log(detail)
                result_lines.append(headline + "\n" + detail)
            if result_lines:
                display = result_lines[0] if len(result_lines) == 1 else f"{len(result_lines)} available results · Open report for Ce³⁺ values"
                self.events.put(("result_summary", display))
            log(text)
            log(f"Results report: {result / 'Results.html'}")
            self.events.put(("result", (result, state, text)))
        except Exception as exc:
            detail = traceback.format_exc()
            log("Run could not complete: " + str(exc))
            log(detail)
            self.events.put(("error", (str(exc), str(log_path))))
        finally:
            self.events.put(("finished", None))

    def _request_stop(self) -> None:
        if not self.active:
            return
        self.stop_event.set()
        self.stop_button.configure(state="disabled")
        self.status.set("Stopping and saving available results; please wait...")
        self.status_label.configure(foreground="#8a6400")
        self._append_log("Stop requested. The sampler will exit at a safe checkpoint and save available results; please wait.")

    def _request_close(self) -> None:
        if self.worker is not None and self.worker.is_alive():
            if self.close_when_finished:
                return
            if messagebox.askyesno("A run is still active", "Stop sampling and exit after available results have been saved?\nChoose No to continue fitting.", parent=self.root):
                self.close_when_finished = True
                self._request_stop()
            return
        self.root.destroy()

    def _append_log(self, text: str) -> None:
        self.log_text.configure(state="normal")
        self.log_text.insert(tk.END, text + "\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state="disabled")

    def _update_progress(self, info: dict) -> None:
        # All displayed counts come from the core, never from timers.
        steps = info.get("steps")
        index = info.get("file_index")
        count = info.get("file_count")
        parts = []
        if index is not None and count is not None:
            parts.append(f"File {index}/{count}")
        if steps is not None:
            parts.append(f"{steps} steps sampled")
        elapsed = info.get("elapsed_seconds")
        if isinstance(elapsed, (int, float)) and math.isfinite(elapsed):
            parts.append(f"Elapsed: {elapsed / 60:.1f} min")
        quantiles = info.get("Ce3_quantiles_percent")
        if isinstance(quantiles, (list, tuple)) and len(quantiles) >= 3:
            try:
                parts.append(f"Ce³⁺ {float(quantiles[1]):.2f}% [{float(quantiles[0]):.2f}, {float(quantiles[2]):.2f}]")
            except (TypeError, ValueError):
                pass
        if parts:
            self.step_status.set(" · ".join(parts))
        total = info.get("total_steps")
        if total is not None and isinstance(steps, (int, float)):
            self.progress_bar.configure(mode="determinate", maximum=total, value=steps)
        elif isinstance(steps, (int, float)) and str(self.progress_bar.cget("mode")) == "determinate":
            self.progress_bar.configure(value=steps)
        if not self.stop_event.is_set():
            self.status.set("Sampling..." if info.get("event") != "complete" else "Preparing results...")

    def _drain_events(self) -> None:
        # Bound each UI pass so a verbose error or batch cannot freeze the window.
        for _ in range(300):
            try:
                kind, payload = self.events.get_nowait()
            except queue.Empty:
                break
            if kind == "log":
                self._append_log(str(payload))
            elif kind == "progress":
                self._update_progress(payload)
            elif kind in {"vamas_inventory", "vamas_error"}:
                self._inspecting_vamas = False
                self.vamas_select_button.configure(state="normal")
                self.run_button.configure(state="normal")
                self.status.set("Ready")
                if kind == "vamas_inventory":
                    path, inventory = payload
                    self._show_vamas_choices(path, inventory)
                else:
                    error, detail = payload
                    self._append_log(detail)
                    messagebox.showerror("Cannot read VAMAS blocks", error, parent=self.root)
            elif kind == "result_summary":
                self.step_status.set(str(payload))
                self.step_status_label.configure(font=("Segoe UI", 11, "bold"), foreground="#153c64")
            elif kind == "result":
                result, state, text = payload
                self.last_output = result
                self.status.set(text)
                color = {"converged": "#18753c", "quick_complete": "#555555", "cancelled": "#8a6400",
                         "failed": "#a32929", "incomplete": "#8a6400"}.get(state, "#555555")
                self.status_label.configure(foreground=color)
            elif kind == "error":
                error, log_path = payload
                self.status.set("Run failed; error details saved in the log.")
                self.status_label.configure(foreground="#a32929")
                if not self.close_when_finished:
                    messagebox.showerror("Fit did not complete", f"{error}\n\nFull error details:\n{log_path}", parent=self.root)
            elif kind == "finished":
                self.progress_bar.stop()
                self._set_running(False)
        if self.close_when_finished and self.worker is not None and not self.worker.is_alive() and self.events.empty():
            self.root.destroy()
            return
        self.root.after(100, self._drain_events)

    def _open_results(self) -> None:
        if self.last_output is not None and self.last_output.is_dir():
            try:
                report = self.last_output / "Results.html"
                os.startfile(str(report if report.is_file() else self.last_output))
            except OSError as exc:
                messagebox.showerror("Cannot open report", str(exc), parent=self.root)


def _all_widgets(widget: tk.Widget):
    for child in widget.winfo_children():
        yield child
        yield from _all_widgets(child)


def _widget_text(widget: tk.Widget) -> str:
    if isinstance(widget, tk.Text):
        return widget.get("1.0", "end-1c")
    try:
        text = str(widget.cget("text"))
        if text:
            return text
    except tk.TclError:
        pass
    try:
        variable = str(widget.cget("textvariable"))
        return str(widget.getvar(variable)) if variable else ""
    except tk.TclError:
        return ""


def _english_interface_snapshot(app: Ce3dXpsFitterApp) -> dict:
    """Inspect application wording while allowing user-provided paths and values."""
    texts = [dict(path="window_title", text=app.root.title())]
    for widget in _all_widgets(app.root):
        # Entries and file lists contain user data, including non-English paths.
        if isinstance(widget, (tk.Entry, ttk.Entry, tk.Listbox, ttk.Combobox)):
            continue
        value = _widget_text(widget)
        if value:
            texts.append(dict(path=str(widget), text=value))
    violations = [item for item in texts if re.search(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]", item["text"])]
    return dict(language="English", initial_log=app.log_text.get("1.0", "end-1c"),
                user_data_excluded="Entry values, selected file paths, and imported spectrum labels may use any language.",
                widget_texts=texts, cjk_violations=violations, passed=not violations)


def _geometry_snapshot(app: Ce3dXpsFitterApp, name: str, size: tuple[int, int], expanded: bool) -> dict:
    """Allocate real Tk geometry with a fully transparent, off-screen window."""
    root = app.root
    if app.advanced_visible != expanded:
        app._toggle_advanced()
    root.geometry(f"{size[0]}x{size[1]}-30000-30000")
    root.deiconify()
    # Processing the initial MapNotify is needed before Tk allocates children.
    root.update()
    root.update_idletasks()
    width, height = root.winfo_width(), root.winfo_height()
    origin_x, origin_y = root.winfo_rootx(), root.winfo_rooty()
    records, hidden = [], []
    for widget in _all_widgets(root):
        if not widget.winfo_ismapped():
            hidden.append(str(widget))
            continue
        x = widget.winfo_rootx() - origin_x
        y = widget.winfo_rooty() - origin_y
        w, h = widget.winfo_width(), widget.winfo_height()
        parent = widget.master
        px = widget.winfo_rootx() - parent.winfo_rootx()
        py = widget.winfo_rooty() - parent.winfo_rooty()
        inside_client = x >= 0 and y >= 0 and x + w <= width and y + h <= height
        inside_parent = (px >= 0 and py >= 0 and px + w <= parent.winfo_width()
                         and py + h <= parent.winfo_height())
        text = _widget_text(widget)
        text_check = dict(checked=False)
        if text and isinstance(widget, (ttk.Label, ttk.Radiobutton, ttk.Button, ttk.LabelFrame)):
            try:
                wraps = int(widget.cget("wraplength") or 0) > 0
            except (tk.TclError, ValueError):
                wraps = False
            if not wraps and "\n" not in text:
                try:
                    font_spec = widget.cget("font")
                except tk.TclError:
                    font_spec = ""
                if not font_spec:
                    style_name = str(widget.cget("style")) or widget.winfo_class()
                    font_spec = ttk.Style(root).lookup(style_name, "font") or "TkDefaultFont"
                font = tkfont.Font(root=root, font=font_spec)
                text_width = font.measure(text)
                # Auto-sized buttons/radios include the theme's actual padding and indicator.
                # Explicit label widths still require independent font measurement.
                required_width = max(text_width, widget.winfo_reqwidth()) if isinstance(widget, (ttk.Button, ttk.Radiobutton)) else text_width + (16 if isinstance(widget, ttk.LabelFrame) else 0)
                text_check = dict(checked=True, measured_width=text_width, required_width=required_width,
                                  fits=required_width <= w)
        records.append(dict(path=str(widget), widget_class=widget.winfo_class(),
                            text=text, x=x, y=y, width=w, height=h,
                            inside_client=inside_client, inside_parent=inside_parent,
                            single_line_text=text_check))
    roles = {}
    for role, widget in (("top_three_rows", app.selection_frame), ("input", app.input_frame),
                         ("output", app.output_frame), ("log", app.log_frame), ("log_text", app.log_text)):
        record = next((item for item in records if item["path"] == str(widget)), None)
        roles[role] = dict(visible=record is not None and record["inside_client"] and record["inside_parent"],
                           height=record["height"] if record else 0)
    result = dict(name=name, requested_client=list(size), actual_client=[width, height],
                  minimum_client=list(root.minsize()), advanced_visible=expanded,
                  requested_content_size=[root.winfo_reqwidth(), root.winfo_reqheight()],
                  widgets=records, intentionally_unmapped_widgets=hidden, essential_areas=roles)
    result["violations"] = [item["path"] for item in records if not item["inside_client"] or not item["inside_parent"]]
    result["text_clipping_violations"] = [item["path"] for item in records
                                           if item["single_line_text"].get("checked") and not item["single_line_text"]["fits"]]
    root.withdraw()
    return result


def run_self_check(report_path: Path, input_path: Path | None = None) -> int:
    """Exercise the frozen runtime without fitting, visible windows, or stdout."""
    report = dict(version=VERSION, status="failed", timestamp=datetime.now().isoformat(),
                  frozen=bool(getattr(sys, "frozen", False)), executable=sys.executable,
                  errors=[], geometry_method="Fully transparent, borderless off-screen mapping only for Tk geometry allocation; no mainloop or sampling.")
    root = None

    def error(stage: str, exc: Exception) -> None:
        report["errors"].append(dict(stage=stage, error=str(exc), traceback=traceback.format_exc()))

    try:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-alpha", 0.0)
        root.overrideredirect(True)
        app = Ce3dXpsFitterApp(root)
        root.update_idletasks()
        report["tk"] = dict(tcl_version=str(root.tk.call("info", "patchlevel")),
                            tk_version=str(root.tk.call("package", "require", "Tk")))
        defaults = dict(model=app.model.get(), satellites=app.satellites.get(), mode=app.mode.get(),
                        walkers=int(app.walkers.get()), steps=int(app.quick_steps.get()), burn=int(app.burn.get()))
        report["defaults"] = defaults
        expected = dict(model="LA", satellites=False, mode="quick", walkers=80, steps=8000, burn=4000)
        if defaults != expected:
            raise AssertionError(f"Unexpected UI defaults: {defaults}")
        top_rows = []
        for row in range(3):
            cells = sorted((widget for widget in app.selection_frame.winfo_children()
                            if int(widget.grid_info().get("row", -1)) == row),
                           key=lambda widget: int(widget.grid_info()["column"]))
            top_rows.append(dict(row=row, label=_widget_text(cells[0]),
                                 options=[dict(text=_widget_text(widget), column=int(widget.grid_info()["column"]),
                                               value=str(widget.cget("value"))) for widget in cells[1:]]))
        report["top_three_rows"] = top_rows
        report["language_check"] = _english_interface_snapshot(app)
        if not report["language_check"]["passed"]:
            raise AssertionError("Non-English application wording was found in the interface.")
        if [item["label"] for item in top_rows] != ["Peak model", "Satellite peaks", "Sampling mode"]:
            raise AssertionError("Top selection row labels do not match the required interface.")
        if any([item["column"] for item in row["options"]] != [1, 2] for row in top_rows):
            raise AssertionError("Each top row must have its two options in columns 1 and 2.")
        lambda_widgets = [str(widget) for widget in _all_widgets(root)
                          if "lambda" in _widget_text(widget).casefold() or "λ" in _widget_text(widget).casefold()]
        report["lambda_ui_absent"] = not lambda_widgets
        if lambda_widgets:
            raise AssertionError(f"Unexpected likelihood weighting control/text: {lambda_widgets}")
        scenes = [
            ("initial_900x800", (900, 800), False),
            ("expanded_900x800", (900, 800), True),
            ("collapsed_minimum", (820, 680), False),
            ("expanded_minimum", (820, 800), True),
        ]
        report["geometry"] = [_geometry_snapshot(app, name, size, expanded) for name, size, expanded in scenes]
        failures = []
        for scene in report["geometry"]:
            if scene["violations"]:
                failures.append(f"{scene['name']}: clipped widgets {scene['violations']}")
            if scene["text_clipping_violations"]:
                failures.append(f"{scene['name']}: clipped text {scene['text_clipping_violations']}")
            for role, state in scene["essential_areas"].items():
                if not state["visible"] or state["height"] < (30 if role == "log_text" else 10):
                    failures.append(f"{scene['name']}: essential area {role} is not usefully visible")
        if failures:
            raise AssertionError("; ".join(failures))
    except Exception as exc:
        error("gui", exc)
    finally:
        if root is not None:
            root.destroy()

    try:
        # Explicit imports also let PyInstaller discover these runtime dependencies.
        import numpy as np
        import scipy
        import pandas as pd
        import emcee
        import matplotlib
        import openpyxl
        import xlrd
        import vamas
        from vamas import Vamas
        import ce3d_xps_fitter_core as core
        from ce3d_xps_model import Ce3dModel
        from ce3d_vamas import inspect_vamas

        if not callable(Vamas) or not callable(core.run_fit):
            raise AssertionError("Required VAMAS/core entry points are unavailable.")
        report["dependencies"] = {name: str(getattr(module, "__version__", "available")) for name, module in (
            ("numpy", np), ("scipy", scipy), ("pandas", pd), ("emcee", emcee), ("matplotlib", matplotlib),
            ("openpyxl", openpyxl), ("xlrd", xlrd), ("vamas", vamas), ("core", core))}
        report["vamas_class_imported"] = True
        energy = np.linspace(870.0, 930.0, 121)
        net = 100.0 + 300.0 * np.exp(-0.5 * ((energy - 900.0) / 5.0) ** 2)
        checks = []
        for shape in ("LA", "Voigt"):
            for satellites in (False, True):
                model = Ce3dModel(energy, net, lineshape=shape, satellites=satellites)
                parameters = model.initial()
                prediction = model.predict(parameters)
                finite = bool(np.isfinite(prediction).all())
                expected_ndim = 15 if satellites else 11
                if model.ndim != expected_ndim or prediction.shape != (1, len(energy)) or not finite:
                    raise AssertionError(f"Invalid {shape} model prediction with satellites={satellites}.")
                log_prob = model.log_prob(parameters)
                if not np.isfinite(log_prob).all():
                    raise AssertionError(f"Non-finite initial posterior: {shape}, satellites={satellites}.")
                checks.append(dict(model=shape, satellites=satellites, ndim=model.ndim,
                                   prediction_shape=list(prediction.shape), prediction_finite=finite,
                                   prediction_min=float(prediction.min()), prediction_max=float(prediction.max()),
                                   initial_log_prob=float(log_prob[0])))
        report["model_checks"] = checks
        if input_path is not None:
            inventory = inspect_vamas(input_path)
            report["vamas_input_check"] = dict(path=str(input_path.resolve()),
                                               block_count=len(inventory["block_inventory"]),
                                               ce3d_candidate_indices=inventory["ce3d_candidate_indices"],
                                               block_inventory=inventory["block_inventory"])
    except Exception as exc:
        error("runtime_dependencies_and_models", exc)
    try:
        from ce3d_export_check import run_export_self_check
        report["export_checks"] = run_export_self_check()
    except Exception as exc:
        error("output_export", exc)
    report["status"] = "failed" if report["errors"] else "passed"
    try:
        report_path = report_path.expanduser().resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    except Exception:
        return 2
    return 0 if report["status"] == "passed" else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Ce 3d XPS Fitter graphical interface")
    parser.add_argument("--self-test", action="store_true", help="Create and destroy a hidden GUI without fitting")
    parser.add_argument("--self-check", type=Path, metavar="JSON_PATH", help="Write hidden GUI and frozen-runtime validation JSON")
    parser.add_argument("--self-check-input", type=Path, metavar="VAMAS_PATH", help="Optionally inspect a VAMAS input during --self-check")
    args = parser.parse_args()
    if args.self_check is not None:
        return run_self_check(args.self_check, args.self_check_input)
    if args.self_test:
        root = tk.Tk()
        root.withdraw()
        app = Ce3dXpsFitterApp(root)
        root.update_idletasks()
        assert app.model.get() == "LA"
        assert app.satellites.get() is False
        assert app.mode.get() == "quick"
        assert app.quick_steps.get() == "8000"
        app.mode.set("convergence")
        app._sync_mode()
        assert str(app.quick_steps_entry.cget("state")) == "disabled"
        app._toggle_advanced()
        root.update_idletasks()
        root.destroy()
        if sys.stdout is not None:
            print("Hidden GUI self-test passed.")
        return 0
    app = Ce3dXpsFitterApp()
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
