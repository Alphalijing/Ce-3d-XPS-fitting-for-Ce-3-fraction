# -*- coding: utf-8 -*-
"""Smart-like local Shirley background used for synthetic benchmarks.

The implementation follows the v6 local Shirley strategy previously found to
closely reproduce the commercial smart background: choose local low-intensity
endpoints from both sides of the Ce 3d window, then run a standard iterative
Shirley calculation without hard clipping the background to the spectrum.
"""
from __future__ import annotations

import numpy as np


def edge_min_endpoints(intensity: np.ndarray, frac: float = 0.10) -> tuple[float, float]:
    """Return endpoint intensities from local minima in both edge windows."""
    y = np.asarray(intensity, dtype=float)
    if y.ndim != 1 or y.size < 12:
        raise ValueError("intensity must be a one-dimensional array with at least 12 points.")
    m = max(6, int(y.size * frac))
    left_idx = int(np.argmin(y[:m]))
    right_idx = int(y.size - m + np.argmin(y[-m:]))
    return float(y[left_idx]), float(y[right_idx])


def shirley_standard(
    binding_energy: np.ndarray,
    intensity: np.ndarray,
    b_start: float,
    b_end: float,
    *,
    max_iter: int = 1400,
    tol: float = 1e-8,
) -> np.ndarray:
    """Compute the standard iterative Shirley background.

    The input may be either ascending or descending in binding energy. No
    max-to-spectrum clipping is applied, so a small amount of crossing is
    allowed, matching the v6 behavior used in the previous smart comparison.
    """
    energy = np.asarray(binding_energy, dtype=float)
    y = np.asarray(intensity, dtype=float)
    if energy.ndim != 1 or y.ndim != 1 or energy.size != y.size:
        raise ValueError("binding_energy and intensity must be one-dimensional arrays of equal length.")
    if energy.size < 2:
        return y.copy()

    rev = False
    if energy[0] > energy[-1]:
        rev = True
        energy = energy[::-1]
        y = y[::-1]
        b_start, b_end = b_end, b_start

    bg = np.linspace(float(b_start), float(b_end), y.size)
    dx = np.diff(energy)
    scale = max(1.0, float(np.nanmax(y) - np.nanmin(y)))

    for _ in range(max_iter):
        signal = y - bg
        trap = 0.5 * (signal[:-1] + signal[1:]) * dx
        cumulative_tail = np.empty(y.size, dtype=float)
        cumulative_tail[-1] = 0.0
        cumulative_tail[:-1] = np.cumsum(trap[::-1])[::-1]
        denom = cumulative_tail[0]

        if abs(denom) < 1e-20:
            candidate = np.linspace(float(b_start), float(b_end), y.size)
        else:
            candidate = float(b_end) + (float(b_start) - float(b_end)) * (cumulative_tail / denom)

        err = float(np.max(np.abs(candidate - bg)))
        bg = 0.7 * candidate + 0.3 * bg
        if err <= tol * scale:
            break

    bg[0] = float(b_start)
    bg[-1] = float(b_end)
    if rev:
        bg = bg[::-1]
    return bg


def smart_like_local_shirley(
    binding_energy: np.ndarray,
    intensity: np.ndarray,
    *,
    edge_frac: float = 0.10,
    max_binding_energy: float | None = 925.0,
) -> tuple[np.ndarray, dict[str, float]]:
    """Generate a smart-like local Shirley background for a Ce 3d spectrum.

    If max_binding_energy is provided, the spectrum is cropped to that upper
    binding-energy boundary before endpoint selection and Shirley iteration,
    matching the v6 workflow for Ce 3d spectra.
    """
    energy = np.asarray(binding_energy, dtype=float)
    y = np.asarray(intensity, dtype=float)
    if max_binding_energy is not None:
        mask = energy <= float(max_binding_energy)
        if int(np.sum(mask)) < 12:
            raise ValueError("Too few data points remain after max_binding_energy cropping.")
        energy = energy[mask]
        y = y[mask]
    b_start, b_end = edge_min_endpoints(y, frac=edge_frac)
    bg = shirley_standard(energy, y, b_start, b_end)
    return bg, {
        "edge_frac": float(edge_frac),
        "max_binding_energy": float(max_binding_energy) if max_binding_energy is not None else float("nan"),
        "b_start": float(b_start),
        "b_end": float(b_end),
        "n_points": float(energy.size),
    }
