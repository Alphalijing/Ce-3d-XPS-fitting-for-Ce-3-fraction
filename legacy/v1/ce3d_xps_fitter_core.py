# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import json
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable

import emcee
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ce3d_xps_model import (
    build_dynamic_confidence_tail_index,
    ce3d_model,
    gaussian_loglike_weighted,
    voigt_profile,
)
from smart_like_shirley_background import smart_like_local_shirley


LogFn = Callable[[str], None]

CE_COLS = [
    "dE",
    "sigma_ce",
    "x",
    "A_ce_total",
    "g3_1",
    "g3_2",
    "g4_1",
    "g4_2",
    "g4_3",
    "log_a_ce",
    "log_b_ce",
]


@dataclass
class PreparedData:
    ce_file: Path
    ce_bg_file: Path
    notes: list[str]


def _log(log: LogFn | None, message: str) -> None:
    if log is not None:
        log(message)


def _safe_label(value: float) -> str:
    return str(value).replace(".", "p").replace("-", "m")


def _read_excel_sheet(path: Path, expected_sheet: str) -> pd.DataFrame:
    xl = pd.ExcelFile(path)
    sheet_map = {str(s).strip().lower(): s for s in xl.sheet_names}
    key = expected_sheet.strip().lower()
    if key not in sheet_map:
        raise ValueError(
            f"{path.name}: sheet '{expected_sheet}' was not found. "
            f"Available sheets: {', '.join(map(str, xl.sheet_names))}"
        )
    return pd.read_excel(path, sheet_name=sheet_map[key], header=None)


def _numeric_scan_from_frame(
    df: pd.DataFrame,
    *,
    sheet_name: str,
    require_background: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, list[str]]:
    notes: list[str] = []
    if df.shape[1] < 3:
        raise ValueError(f"{sheet_name}: sheet must contain at least columns A and C.")
    if require_background and df.shape[1] < 5:
        raise ValueError(f"{sheet_name}: provided-background mode requires columns A, C, and E.")

    energy = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    intensity = pd.to_numeric(df.iloc[:, 2], errors="coerce")
    if require_background:
        background = pd.to_numeric(df.iloc[:, 4], errors="coerce")
        mask = energy.notna() & intensity.notna() & background.notna()
    else:
        background = None
        mask = energy.notna() & intensity.notna()
    if int(mask.sum()) < 20:
        raise ValueError(f"{sheet_name}: fewer than 20 numeric data rows were detected.")

    e = energy[mask].to_numpy(dtype=float)
    y = intensity[mask].to_numpy(dtype=float)
    bg = background[mask].to_numpy(dtype=float) if background is not None else None

    if bg is not None and bg.size >= 2 and np.isclose(bg[0], 0.0) and not np.isclose(bg[1], 0.0):
        e = e[1:]
        y = y[1:]
        bg = bg[1:]
        notes.append(f"{sheet_name}: dropped the first data row because the background value is 0.")

    order = np.argsort(e)
    return e[order], y[order], bg[order] if bg is not None else None, notes


def _read_processed_file(path: Path, *, require_background: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    suffix = path.suffix.lower()
    if suffix in {".xls", ".xlsx"}:
        df = pd.read_excel(path, header=None)
    elif suffix in {".csv", ".txt", ".dat"}:
        df = pd.read_csv(path, header=None)
    else:
        raise ValueError(f"{path.name}: unsupported file type. Use .csv, .txt, .dat, .xls, or .xlsx.")

    if require_background and df.shape[1] < 3:
        raise ValueError(f"{path.name}: processed input must have three columns: energy, intensity, background.")
    if not require_background and df.shape[1] < 2:
        raise ValueError(f"{path.name}: processed input must have at least two columns: energy, intensity.")

    energy = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    intensity = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    if require_background:
        background = pd.to_numeric(df.iloc[:, 2], errors="coerce")
        mask = energy.notna() & intensity.notna() & background.notna()
    else:
        background = None
        mask = energy.notna() & intensity.notna()
    if int(mask.sum()) < 20:
        raise ValueError(f"{path.name}: fewer than 20 numeric data rows were detected.")

    e = energy[mask].to_numpy(dtype=float)
    y = intensity[mask].to_numpy(dtype=float)
    bg = background[mask].to_numpy(dtype=float) if background is not None else None
    order = np.argsort(e)
    return e[order], y[order], bg[order] if bg is not None else None


def _looks_like_ce3d(path: Path, e: np.ndarray) -> bool:
    med = float(np.nanmedian(e))
    name = path.name.lower()
    return 850.0 <= med <= 940.0 or "ce" in name or "3d" in name


def _write_two_col(path: Path, e: np.ndarray, y: np.ndarray) -> None:
    pd.DataFrame({"Binding Energy (eV)": e, "Intensity (a.u.)": y}).to_csv(
        path, index=False, header=False, encoding="utf-8-sig"
    )


def _write_three_col(path: Path, e: np.ndarray, y: np.ndarray, bg: np.ndarray) -> None:
    pd.DataFrame(
        {
            "Binding Energy (eV)": e,
            "Intensity (a.u.)": y,
            "Background (a.u.)": bg,
        }
    ).to_csv(path, index=False, encoding="utf-8-sig")


def _compute_v6_background(e: np.ndarray, y_raw: np.ndarray, *, log: LogFn | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Compute the v6 smart-like local Shirley background and return cropped arrays."""
    e = np.asarray(e, dtype=float)
    y_raw = np.asarray(y_raw, dtype=float)
    mask = e <= 925.0
    if int(np.sum(mask)) >= 20:
        e_use = e[mask]
        y_use = y_raw[mask]
    else:
        e_use = e
        y_use = y_raw
    bg, info = smart_like_local_shirley(e_use, y_use, edge_frac=0.10, max_binding_energy=None)
    _log(log, "No Shirley/smart background was provided; generated v6 smart-like local Shirley background.")
    return e_use, y_use, bg, info


def prepare_input_files(
    input_paths: Iterable[str | Path],
    *,
    raw_xps: bool,
    background_provided: bool,
    outdir: str | Path,
    log: LogFn | None = None,
) -> PreparedData:
    paths = [Path(p) for p in input_paths]
    if not paths:
        raise ValueError("No input file was provided.")
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(str(p))

    outdir = Path(outdir)
    prepared = outdir / "prepared_input"
    prepared.mkdir(parents=True, exist_ok=True)
    notes: list[str] = []

    if raw_xps:
        if len(paths) != 1:
            raise ValueError("Raw XPS mode expects one Excel workbook.")
        book = paths[0]
        if book.suffix.lower() not in {".xls", ".xlsx"}:
            raise ValueError("Raw XPS mode requires .xls or .xlsx input.")
        _log(log, f"Reading raw XPS workbook: {book}")
        _read_excel_sheet(book, "XPS Survey")
        ce_df = _read_excel_sheet(book, "Ce3d Scan")
        e_ce, y_ce_raw, bg_ce, ce_notes = _numeric_scan_from_frame(
            ce_df,
            sheet_name="Ce3d Scan",
            require_background=background_provided,
        )
        notes.extend(ce_notes)
        if bg_ce is None:
            e_ce, y_ce_raw, bg_ce, bg_info = _compute_v6_background(e_ce, y_ce_raw, log=log)
            notes.append("Ce3d Scan: background generated by v6 smart-like local Shirley.")
            notes.append(f"v6 background info: {bg_info}")
        ce_file = prepared / "processed_ce3d_bg_sub.csv"
        ce_bg_file = prepared / "ce3d_background.csv"
        _write_two_col(ce_file, e_ce, y_ce_raw - bg_ce)
        _write_two_col(ce_bg_file, e_ce, bg_ce)
        _write_three_col(prepared / "ce3d_raw_with_background.csv", e_ce, y_ce_raw, bg_ce)
        shutil.copy2(book, prepared / book.name)
        return PreparedData(ce_file=ce_file, ce_bg_file=ce_bg_file, notes=notes)

    ce_file = None
    ce_bg_file = None
    for path in paths:
        _log(log, f"Reading processed data file: {path}")
        e, y_raw, bg = _read_processed_file(path, require_background=background_provided)
        if not _looks_like_ce3d(path, e):
            _log(log, f"Skipping non-Ce 3d file: {path.name}")
            continue
        if bg is None:
            e, y_raw, bg, bg_info = _compute_v6_background(e, y_raw, log=log)
            notes.append(f"{path.name}: background generated by v6 smart-like local Shirley.")
            notes.append(f"v6 background info: {bg_info}")
        ce_file = prepared / "processed_ce3d_bg_sub.csv"
        ce_bg_file = prepared / "ce3d_background.csv"
        _write_two_col(ce_file, e, y_raw - bg)
        _write_two_col(ce_bg_file, e, bg)
        _write_three_col(prepared / "ce3d_raw_with_background.csv", e, y_raw, bg)
        break

    if ce_file is None or ce_bg_file is None:
        raise ValueError("Processed mode requires a Ce 3d file.")
    return PreparedData(ce_file=ce_file, ce_bg_file=ce_bg_file, notes=notes)


def read_spectrum_2col(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, header=None)
    if df.shape[1] < 2:
        raise ValueError(f"{path}: expected two columns.")
    e = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    y = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    mask = e.notna() & y.notna()
    e_arr = e[mask].to_numpy(dtype=float)
    y_arr = y[mask].to_numpy(dtype=float)
    order = np.argsort(e_arr)
    return e_arr[order], y_arr[order]


def _ce_param_dict(theta: np.ndarray) -> dict[str, float]:
    return {
        "dE": theta[0],
        "sigma_ce": theta[1],
        "x": theta[2],
        "A_ce_total": theta[3],
        "g3_1": theta[4],
        "g3_2": theta[5],
        "g4_1": theta[6],
        "g4_2": theta[7],
        "g4_3": theta[8],
        "log_a_ce": theta[9],
        "log_b_ce": theta[10],
    }


def _metrics(y: np.ndarray, yhat: np.ndarray) -> dict[str, float]:
    resid = y - yhat
    sse = float(np.sum(resid**2))
    rmse = float(np.sqrt(np.mean(resid**2)))
    sst = float(np.sum((y - np.mean(y)) ** 2)) + 1e-15
    return {"RMSE": rmse, "R2": float(1.0 - sse / sst), "SSE": sse}


def _ce_log_prior(theta: np.ndarray, I_ce: np.ndarray) -> float:
    p = _ce_param_dict(theta)
    dE = float(p["dE"])
    sigma_ce = float(p["sigma_ce"])
    x = float(p["x"])
    A_total = float(p["A_ce_total"])
    gammas = [float(p[k]) for k in ["g3_1", "g3_2", "g4_1", "g4_2", "g4_3"]]
    log_a_ce = float(p["log_a_ce"])
    log_b_ce = float(p["log_b_ce"])

    if not (-1.0 < dE < 1.0):
        return -np.inf
    if not (0.05 < sigma_ce < 1.5):
        return -np.inf
    if not (0.0 < x < 1.0):
        return -np.inf
    amp_max_ce = max(float(np.max(I_ce)) * 200.0, 1e3)
    if not (1.0 < A_total < amp_max_ce):
        return -np.inf
    if any(not (0.10 < g < 2.0) for g in gammas):
        return -np.inf
    if not (-20.0 < log_a_ce < 20.0 and -20.0 < log_b_ce < 20.0):
        return -np.inf

    lp = 0.0
    lp += -0.5 * ((sigma_ce - 0.5) / 0.2) ** 2
    lp += -0.5 * ((dE - 0.0) / 0.2) ** 2
    for g in gammas:
        lp += -0.5 * ((g - 0.5) / 0.3) ** 2
    lp += -0.5 * ((log_a_ce - np.log(max(np.var(I_ce) * 0.05, 1e-6))) / 2.0) ** 2
    lp += -0.5 * ((log_b_ce - np.log(1e-2)) / 2.0) ** 2
    return float(lp)


def _ce_log_posterior(theta: np.ndarray, E_ce: np.ndarray, I_ce: np.ndarray, C_ce: np.ndarray) -> float:
    lp = _ce_log_prior(theta, I_ce)
    if not np.isfinite(lp):
        return -np.inf
    p = _ce_param_dict(theta)
    y_ce, _ = ce3d_model(E_ce, p)
    if np.any(y_ce <= 0):
        return -np.inf
    ll = gaussian_loglike_weighted(I_ce, y_ce, float(p["log_a_ce"]), float(p["log_b_ce"]), C_ce)
    return float(lp + ll) if np.isfinite(ll) else -np.inf


def _ce_initial(E_ce: np.ndarray, I_ce: np.ndarray) -> np.ndarray:
    ce_area = max(float(np.trapz(np.clip(I_ce, 0, None), E_ce)), 100.0)
    return np.array(
        [
            0.0,
            0.55,
            0.2,
            ce_area,
            0.45,
            0.55,
            0.45,
            0.55,
            0.65,
            np.log(max(np.var(I_ce) * 0.05, 1e-6)),
            np.log(1e-2),
        ],
        dtype=float,
    )


def _run_sampler(
    initial: np.ndarray,
    args: tuple,
    *,
    walkers: int,
    steps: int,
    burn: int,
    thin: int,
    seed: int,
) -> tuple[emcee.EnsembleSampler, np.ndarray]:
    ndim = initial.size
    if walkers < 2 * ndim:
        raise ValueError(f"walkers must be >= {2 * ndim} for ndim={ndim}.")
    rng = np.random.default_rng(seed)
    jitter = np.where(np.abs(initial) > 1e-8, 0.03 * np.abs(initial), 0.02)
    p0 = initial + jitter * rng.standard_normal((walkers, ndim))
    sampler = emcee.EnsembleSampler(walkers, ndim, _ce_log_posterior, args=args)
    sampler.run_mcmc(p0, steps, progress=False)
    flat = sampler.get_chain(discard=burn, thin=thin, flat=True)
    return sampler, flat


def _param_summary(flat: np.ndarray) -> pd.DataFrame:
    rows = []
    for i, col in enumerate(CE_COLS):
        rows.append(
            {
                "parameter": col,
                "median": float(np.median(flat[:, i])),
                "mean": float(np.mean(flat[:, i])),
                "std": float(np.std(flat[:, i])),
                "q025": float(np.quantile(flat[:, i], 0.025)),
                "q975": float(np.quantile(flat[:, i], 0.975)),
            }
        )
    return pd.DataFrame(rows)


def _plot_ce(outdir: Path, E_ce: np.ndarray, I_ce: np.ndarray, y_ce: np.ndarray, comps: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(E_ce, I_ce, color="black", lw=1.0, label="Ce 3d data")
    ax.plot(E_ce, y_ce, color="#235789", lw=1.5, label="Total fit")
    for comp in comps:
        color = "#53B0A4" if comp["component"].startswith("Ce3") else "#EA7E63"
        curve = comp["area"] * voigt_profile(E_ce, comp["center_eV"], comp["sigma"], comp["gamma"])
        ax.fill_between(E_ce, curve, 0, color=color, alpha=0.20, lw=0)
        ax.plot(E_ce, curve, color=color, lw=0.8, alpha=0.9)
    ax.set_xlabel("Binding Energy (eV)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.invert_xaxis()
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "fit_ce3d.png", dpi=180)
    plt.close(fig)


def fit_ce3d(
    ce_file: Path,
    ce_bg_file: Path,
    outdir: Path,
    *,
    lambda_ce: float,
    walkers: int,
    steps: int,
    burn: int,
    thin: int,
    seed: int,
    log: LogFn | None = None,
) -> dict:
    _log(log, f"Running Ce 3d MCMC, lambda={lambda_ce:g}")
    E_ce, I_ce = read_spectrum_2col(ce_file)
    E_bg, I_bg = read_spectrum_2col(ce_bg_file)
    tail_index, tail_info = build_dynamic_confidence_tail_index(E_ce, I_ce, E_bg, I_bg)
    confidence = 1.0 / (1.0 + float(lambda_ce) * tail_index)
    _log(log, f"Dynamic confidence reset energies: {tail_info.get('dynamic_reset_energies_eV', '')}")
    sampler, flat = _run_sampler(
        _ce_initial(E_ce, I_ce),
        (E_ce, I_ce, confidence),
        walkers=walkers,
        steps=steps,
        burn=burn,
        thin=thin,
        seed=seed,
    )
    pd.DataFrame(flat, columns=CE_COLS).to_csv(outdir / "posterior_samples_ce3d.csv", index=False, encoding="utf-8-sig")
    theta = np.median(flat, axis=0)
    p = _ce_param_dict(theta)
    y_ce, comps = ce3d_model(E_ce, p)
    metrics = _metrics(I_ce, y_ce)
    pd.DataFrame(
        {
            "E_ce": E_ce,
            "I_ce": I_ce,
            "I_ce_fit": y_ce,
            "resid_ce": I_ce - y_ce,
            "tail_index_ce": tail_index,
            "confidence_ce": confidence,
        }
    ).to_csv(outdir / "fit_curve_ce3d.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(comps).to_csv(outdir / "peak_table_ce3d.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(comps).to_csv(outdir / "peak_table.csv", index=False, encoding="utf-8-sig")
    _param_summary(flat).to_csv(outdir / "posterior_param_summary_ce3d.csv", index=False, encoding="utf-8-sig")
    _plot_ce(outdir, E_ce, I_ce, y_ce, comps)
    return {
        "lambda_ce": float(lambda_ce),
        "Ce3_fraction_mean": float(np.mean(flat[:, 2])),
        "Ce3_fraction_std": float(np.std(flat[:, 2])),
        "Ce3_fraction_median": float(np.median(flat[:, 2])),
        "acceptance_fraction_mean_ce3d": float(np.mean(sampler.acceptance_fraction)),
        "Ce3d_metrics": metrics,
        "posterior_median_ce3d": {col: float(val) for col, val in zip(CE_COLS, theta)},
        "tail_info": tail_info,
    }


def run_fit(
    input_paths: Iterable[str | Path],
    *,
    raw_xps: bool,
    lambda_ce: float,
    background_provided: bool = True,
    outdir: str | Path,
    walkers: int = 80,
    steps: int = 5000,
    burn: int = 1500,
    thin: int = 10,
    seed: int = 42,
    log: LogFn | None = None,
) -> Path:
    outroot = Path(outdir)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = outroot / f"ce3d_fit_lambda_{_safe_label(float(lambda_ce))}_{stamp}"
    run_dir.mkdir(parents=True, exist_ok=False)

    prepared = prepare_input_files(
        input_paths,
        raw_xps=raw_xps,
        background_provided=background_provided,
        outdir=run_dir,
        log=log,
    )
    for note in prepared.notes:
        _log(log, note)

    summary = {
        "raw_xps_input": bool(raw_xps),
        "background_provided": bool(background_provided),
        "fit_target": "Ce 3d only",
        "lambda_ce": float(lambda_ce),
        "input_files": [str(Path(p)) for p in input_paths],
        "prepared_ce_file": str(prepared.ce_file),
        "prepared_ce_bg_file": str(prepared.ce_bg_file),
        "mcmc": {"walkers": walkers, "steps": steps, "burn": burn, "thin": thin, "seed": seed},
        "notes": prepared.notes,
    }
    summary.update(
        fit_ce3d(
            prepared.ce_file,
            prepared.ce_bg_file,
            run_dir,
            lambda_ce=lambda_ce,
            walkers=walkers,
            steps=steps,
            burn=burn,
            thin=thin,
            seed=seed,
            log=log,
        )
    )

    with (run_dir / "fit_result.json").open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    _log(log, f"Done. Output directory: {run_dir}")
    return run_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Automated physically constrained Ce 3d XPS fitting.")
    parser.add_argument("--input", nargs="+", required=True, help="Input raw XPS workbook or processed Ce 3d files.")
    parser.add_argument("--raw-xps", action="store_true", help="Input is a raw XPS workbook.")
    parser.add_argument(
        "--auto-shirley",
        action="store_true",
        help="No Shirley/smart background is provided; compute v6 smart-like local Shirley background.",
    )
    parser.add_argument("--lambda-ce", type=float, default=2.0)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--walkers", type=int, default=80)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--burn", type=int, default=1500)
    parser.add_argument("--thin", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_fit(
        args.input,
        raw_xps=args.raw_xps,
        background_provided=not args.auto_shirley,
        lambda_ce=args.lambda_ce,
        outdir=args.outdir,
        walkers=args.walkers,
        steps=args.steps,
        burn=args.burn,
        thin=args.thin,
        seed=args.seed,
        log=print,
    )


if __name__ == "__main__":
    main()
