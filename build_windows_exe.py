# -*- coding: utf-8 -*-
from __future__ import annotations

import importlib._bootstrap_external as bootstrap_external
import os
import sys
from datetime import datetime


def direct_write_pyc(path, data, mode=0o666):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)
    try:
        os.chmod(path, mode)
    except OSError:
        pass


bootstrap_external._write_atomic = direct_write_pyc

stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

sys.argv = [
    "pyinstaller",
    "--noconfirm",
    "--clean",
    "--windowed",
    "--noupx",
    "--collect-data",
    "tkinterdnd2",
    "--name",
    "Ce3d_XPS_Fitter",
    "--workpath",
    f"build_patched_{stamp}",
    "--distpath",
    "dist_patched",
    "ce3d_xps_fitter_gui.py",
]

from PyInstaller.__main__ import run

run()
