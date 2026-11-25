# -*- coding: utf-8 -*-
from __future__ import annotations

import queue
import threading
import traceback
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
except Exception:
    TkinterDnD = None
    DND_FILES = None

from ce3d_xps_fitter_core import run_fit


class Ce3dXpsFitterApp:
    def __init__(self) -> None:
        self.dnd_available = False
        if TkinterDnD is not None:
            try:
                self.root = TkinterDnD.Tk()
                self.dnd_available = True
            except Exception:
                self.root = tk.Tk()
        else:
            self.root = tk.Tk()

        self.root.title("Ce 3d XPS Fitter")
        self.root.geometry("780x650")
        self.root.minsize(720, 610)

        self.raw_xps = tk.IntVar(value=1)
        self.background_provided = tk.IntVar(value=1)
        self.lambda_choice = tk.IntVar(value=2)
        self.lambda_custom = tk.StringVar(value="2")
        self.output_dir = tk.StringVar(value=str(Path.cwd() / "ce3d_xps_fit_outputs"))
        self.input_files: list[Path] = []
        self.log_queue: queue.Queue[str] = queue.Queue()
        self.worker: threading.Thread | None = None

        self._build_ui()
        self.root.after(150, self._drain_log_queue)

    def mainloop(self) -> None:
        self.root.mainloop()

    def _build_ui(self) -> None:
        frame = ttk.Frame(self.root, padding=14)
        frame.pack(fill=tk.BOTH, expand=True)
        frame.columnconfigure(0, weight=1)

        title = ttk.Label(frame, text="Automated Ce 3d XPS Fitting", font=("Segoe UI", 14, "bold"))
        title.grid(row=0, column=0, sticky="w", pady=(0, 10))

        input_box = ttk.LabelFrame(frame, text="1. Input data type", padding=10)
        input_box.grid(row=1, column=0, sticky="ew", pady=5)
        ttk.Radiobutton(
            input_box,
            text="Raw XPS workbook (.xls/.xlsx)",
            variable=self.raw_xps,
            value=1,
        ).grid(row=0, column=0, sticky="w", padx=4, pady=2)
        ttk.Radiobutton(
            input_box,
            text="Processed Ce 3d file (CSV/table; see README)",
            variable=self.raw_xps,
            value=0,
        ).grid(row=1, column=0, sticky="w", padx=4, pady=2)

        bg_box = ttk.LabelFrame(frame, text="2. Shirley/smart background", padding=10)
        bg_box.grid(row=2, column=0, sticky="ew", pady=5)
        ttk.Radiobutton(
            bg_box,
            text="Provided in input file",
            variable=self.background_provided,
            value=1,
        ).grid(row=0, column=0, sticky="w", padx=4, pady=2)
        ttk.Radiobutton(
            bg_box,
            text="Not provided; compute v6 smart-like Shirley background",
            variable=self.background_provided,
            value=0,
        ).grid(row=1, column=0, sticky="w", padx=4, pady=2)

        lambda_box = ttk.LabelFrame(frame, text="3. Ce 3d likelihood weighting coefficient", padding=10)
        lambda_box.grid(row=3, column=0, sticky="ew", pady=5)
        lambda_box.columnconfigure(3, weight=1)
        ttk.Radiobutton(
            lambda_box,
            text="Unweighted (lambda = 0)",
            variable=self.lambda_choice,
            value=0,
            command=self._sync_lambda_entry,
        ).grid(row=0, column=0, sticky="w", padx=4, pady=2)
        ttk.Radiobutton(
            lambda_box,
            text="lambda = 2 (recommended)",
            variable=self.lambda_choice,
            value=2,
            command=self._sync_lambda_entry,
        ).grid(row=0, column=1, sticky="w", padx=18, pady=2)
        ttk.Radiobutton(
            lambda_box,
            text="Custom lambda:",
            variable=self.lambda_choice,
            value=3,
            command=self._sync_lambda_entry,
        ).grid(row=0, column=2, sticky="w", padx=(18, 2), pady=2)
        validate = (self.root.register(self._validate_float_text), "%P")
        self.lambda_entry = ttk.Entry(lambda_box, textvariable=self.lambda_custom, width=12, validate="key", validatecommand=validate)
        self.lambda_entry.grid(row=0, column=3, sticky="w", padx=2, pady=2)
        self._sync_lambda_entry()

        file_box = ttk.LabelFrame(frame, text="4. Input file", padding=10)
        file_box.grid(row=4, column=0, sticky="nsew", pady=5)
        file_box.columnconfigure(0, weight=1)
        file_box.rowconfigure(1, weight=1)
        ttk.Label(file_box, text="Drop the XPS data file here, or use the Browse button.").grid(
            row=0, column=0, sticky="w", pady=(0, 5)
        )
        self.file_list = tk.Listbox(file_box, height=7, activestyle="none")
        self.file_list.grid(row=1, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(file_box, orient=tk.VERTICAL, command=self.file_list.yview)
        scroll.grid(row=1, column=1, sticky="ns")
        self.file_list.configure(yscrollcommand=scroll.set)
        if self.dnd_available:
            self.file_list.drop_target_register(DND_FILES)
            self.file_list.dnd_bind("<<Drop>>", self._on_drop)

        button_row = ttk.Frame(file_box)
        button_row.grid(row=2, column=0, sticky="w", pady=(8, 0))
        ttk.Button(button_row, text="Browse file...", command=self._browse_files).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(button_row, text="Clear", command=self._clear_files).pack(side=tk.LEFT)

        out_box = ttk.LabelFrame(frame, text="5. Output folder", padding=10)
        out_box.grid(row=5, column=0, sticky="ew", pady=5)
        out_box.columnconfigure(0, weight=1)
        ttk.Entry(out_box, textvariable=self.output_dir).grid(row=0, column=0, sticky="ew", padx=(0, 8))
        ttk.Button(out_box, text="Browse...", command=self._browse_output).grid(row=0, column=1)

        run_row = ttk.Frame(frame)
        run_row.grid(row=6, column=0, sticky="ew", pady=(8, 6))
        self.run_button = ttk.Button(run_row, text="Start fitting", command=self._start_fit)
        self.run_button.pack(side=tk.LEFT)
        ttk.Label(run_row, text="Default MCMC: walkers=80, steps=5000, burn=1500, thin=10").pack(side=tk.LEFT, padx=14)

        log_box = ttk.LabelFrame(frame, text="6. Run log / error messages", padding=10)
        log_box.grid(row=7, column=0, sticky="nsew", pady=5)
        frame.rowconfigure(7, weight=1)
        log_box.columnconfigure(0, weight=1)
        log_box.rowconfigure(0, weight=1)
        self.log_text = tk.Text(log_box, height=10, wrap="word")
        self.log_text.grid(row=0, column=0, sticky="nsew")
        log_scroll = ttk.Scrollbar(log_box, orient=tk.VERTICAL, command=self.log_text.yview)
        log_scroll.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=log_scroll.set)
        if not self.dnd_available:
            self._append_log("Drag-and-drop support is unavailable on this system. Use Browse file instead.")

    def _validate_float_text(self, text: str) -> bool:
        if text.strip() == "":
            return True
        try:
            value = float(text)
        except ValueError:
            return False
        return value >= 0

    def _sync_lambda_entry(self) -> None:
        choice = self.lambda_choice.get()
        if choice == 0:
            self.lambda_custom.set("0")
            self.lambda_entry.configure(state="disabled")
        elif choice == 2:
            self.lambda_custom.set("2")
            self.lambda_entry.configure(state="disabled")
        else:
            self.lambda_entry.configure(state="normal")
            self.lambda_entry.focus_set()

    def _on_drop(self, event) -> None:
        files = self.root.tk.splitlist(event.data)
        self._add_files(files)

    def _browse_files(self) -> None:
        filenames = filedialog.askopenfilenames(
            title="Select XPS data file",
            filetypes=[("XPS/CSV files", "*.xls *.xlsx *.csv *.txt *.dat"), ("All files", "*.*")],
        )
        self._add_files(filenames)

    def _add_files(self, files) -> None:
        for item in files:
            path = Path(str(item).strip("{}"))
            if path.exists() and path not in self.input_files:
                self.input_files.append(path)
                self.file_list.insert(tk.END, str(path))

    def _clear_files(self) -> None:
        self.input_files.clear()
        self.file_list.delete(0, tk.END)

    def _browse_output(self) -> None:
        dirname = filedialog.askdirectory(title="Select output folder", initialdir=self.output_dir.get() or str(Path.cwd()))
        if dirname:
            self.output_dir.set(dirname)

    def _get_lambda(self) -> float:
        if self.lambda_choice.get() == 0:
            return 0.0
        if self.lambda_choice.get() == 2:
            return 2.0
        text = self.lambda_custom.get().strip()
        if not text:
            raise ValueError("Custom lambda cannot be empty.")
        value = float(text)
        if value < 0:
            raise ValueError("lambda must be a non-negative number.")
        return value

    def _start_fit(self) -> None:
        if self.worker is not None and self.worker.is_alive():
            messagebox.showinfo("Running", "A fitting task is already running.")
            return
        try:
            if not self.input_files:
                raise ValueError("Please select or drop an input file first.")
            lambda_ce = self._get_lambda()
            outdir = Path(self.output_dir.get())
            outdir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            messagebox.showerror("Input error", str(exc))
            return

        self.run_button.configure(state="disabled")
        self._append_log("")
        self._append_log("========== New Ce 3d fit ==========")
        self.worker = threading.Thread(
            target=self._fit_worker,
            args=(
                list(self.input_files),
                bool(self.raw_xps.get()),
                bool(self.background_provided.get()),
                lambda_ce,
                outdir,
            ),
            daemon=True,
        )
        self.worker.start()

    def _fit_worker(
        self,
        input_files: list[Path],
        raw_xps: bool,
        background_provided: bool,
        lambda_ce: float,
        outdir: Path,
    ) -> None:
        try:
            run_dir = run_fit(
                input_files,
                raw_xps=raw_xps,
                background_provided=background_provided,
                lambda_ce=lambda_ce,
                outdir=outdir,
                log=lambda msg: self.log_queue.put(str(msg)),
            )
            self.log_queue.put(f"Output completed: {run_dir}")
            self.log_queue.put("__DONE__")
        except Exception:
            self.log_queue.put(traceback.format_exc())
            self.log_queue.put("__FAILED__")

    def _append_log(self, message: str) -> None:
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)

    def _drain_log_queue(self) -> None:
        try:
            while True:
                msg = self.log_queue.get_nowait()
                if msg == "__DONE__":
                    self.run_button.configure(state="normal")
                    messagebox.showinfo("Done", "Ce 3d fitting completed.")
                elif msg == "__FAILED__":
                    self.run_button.configure(state="normal")
                    messagebox.showerror("Fitting failed", "The fitting task failed. See the log box.")
                else:
                    self._append_log(msg)
        except queue.Empty:
            pass
        self.root.after(150, self._drain_log_queue)


def main() -> None:
    app = Ce3dXpsFitterApp()
    app.mainloop()


if __name__ == "__main__":
    main()
