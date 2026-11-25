# -*- coding: utf-8 -*-
"""Build a windowed one-directory Windows distribution using this interpreter."""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


APP_NAME = "Ce3d_XPS_Fitter_v2.0.1_EN"


def main() -> int:
    base = Path(__file__).resolve().parent
    entry = base / "ce3d_xps_fitter_gui.py"
    if not entry.is_file() or not (base / "ce3d_xps_fitter_core.py").is_file():
        print("Build files are incomplete: GUI and core must be next to this script.", file=sys.stderr)
        return 2
    if importlib.util.find_spec("PyInstaller") is None:
        print(f'PyInstaller is missing. Install it with: "{sys.executable}" -m pip install pyinstaller', file=sys.stderr)
        return 2
    command = [
        sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--onedir", "--windowed", "--noupx",
        "--name", APP_NAME, "--distpath", str(base / "dist"),
        "--workpath", str(base / "build" / APP_NAME), "--specpath", str(base / "build"),
        "--paths", str(base), "--hidden-import", "matplotlib.backends.backend_agg",
        "--hidden-import", "matplotlib.backends.backend_svg",
        "--hidden-import", "xlrd", "--hidden-import", "openpyxl", "--hidden-import", "vamas",
    ]
    for unused in ("IPython", "pytest", "notebook", "PyQt5", "PyQt6", "PySide2", "PySide6", "torch", "tensorflow"):
        command.extend(["--exclude-module", unused])
    if importlib.util.find_spec("tkinterdnd2") is not None:
        command.extend(["--collect-data", "tkinterdnd2"])
        print("Including optional drag-and-drop support.")
    else:
        print("tkinterdnd2 is not installed; the application will use the file browser.")
    command.append(str(entry))
    print(f"Building {APP_NAME} with {sys.executable}")
    completed = subprocess.run(command, cwd=base, check=False)
    if completed.returncode == 0:
        print(f"Ready: {base / 'dist' / APP_NAME / (APP_NAME + '.exe')}")
        print("Distribute the ENTIRE application folder, including _internal; the exe alone is insufficient.")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
