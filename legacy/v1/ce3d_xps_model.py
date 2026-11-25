# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
from scipy.signal import savgol_filter
from scipy.special import wofz


DELTA_SO = 18.6

# Ce 3d 5/2 base positions in eV.
E4_BASE = np.array([882.6, 888.8, 898.4], dtype=float)
E3_BASE = np.array([880.6, 885.2], dtype=float)

# Pair-level multiplet templates.
W4_PAIR = np.array([3.9, 2.5, 3.8], dtype=float)
W3_PAIR = np.array([2.9, 5.4], dtype=float)
W4_PAIR = W4_PAIR / W4_PAIR.sum()
W3_PAIR = W3_PAIR / W3_PAIR.sum()

# Degeneracy-derived spin-orbit intensity ratio:
# 3d5/2 : 3d3/2 = 3 : 2.
R_52 = 3.0 / 5.0
R_32 = 2.0 / 5.0


def voigt_profile(x: np.ndarray, mu: float, sigma: float, gamma: float) -> np.ndarray:
    """Unit-area Voigt profile."""
    z = ((x - mu) + 1j * gamma) / (sigma * np.sqrt(2.0))
    return np.real(wofz(z)) / (sigma * np.sqrt(2.0 * np.pi))


def voigt_fwhm_approx(sigma: float, gamma: float) -> float:
    """Approximate Voigt full width at half maximum."""
    f_l = 2.0 * gamma
    f_g = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma
    return float(0.5346 * f_l + np.sqrt(0.2166 * f_l * f_l + f_g * f_g))


def ce3d_model(E_ce: np.ndarray, p: Dict[str, np.ndarray | float]) -> Tuple[np.ndarray, List[Dict[str, float]]]:
    """
    Physically constrained ten-component Ce 3d model.

    The model uses:
    - fixed 18.6 eV spin-orbit splitting,
    - fixed 3d5/2:3d3/2 degeneracy ratio of 3:2,
    - fixed Ce3+/Ce4+ multiplet templates,
    - one shared Gaussian broadening,
    - grouped Lorentzian broadenings.
    """
    dE = float(p["dE"])
    sigma = float(p["sigma_ce"])
    x = float(p["x"])
    A_total = float(p["A_ce_total"])

    A3 = x * A_total
    A4 = (1.0 - x) * A_total

    g3 = [float(p["g3_1"]), float(p["g3_2"])]
    g4 = [float(p["g4_1"]), float(p["g4_2"]), float(p["g4_3"])]

    y = np.zeros_like(E_ce, dtype=float)
    comps: list[dict[str, float]] = []

    for i, (e5, w, g) in enumerate(zip(E3_BASE, W3_PAIR, g3), start=1):
        c52 = (A3 * w * R_52) * voigt_profile(E_ce, e5 + dE, sigma, g)
        c32 = (A3 * w * R_32) * voigt_profile(E_ce, e5 + DELTA_SO + dE, sigma, g)
        y += c52 + c32
        comps.append(
            {
                "spectrum": "Ce3d",
                "component": f"Ce3_pair{i}_3d5/2",
                "center_eV": float(e5 + dE),
                "sigma": sigma,
                "gamma": g,
                "fwhm_approx": voigt_fwhm_approx(sigma, g),
                "area": float(A3 * w * R_52),
            }
        )
        comps.append(
            {
                "spectrum": "Ce3d",
                "component": f"Ce3_pair{i}_3d3/2",
                "center_eV": float(e5 + DELTA_SO + dE),
                "sigma": sigma,
                "gamma": g,
                "fwhm_approx": voigt_fwhm_approx(sigma, g),
                "area": float(A3 * w * R_32),
            }
        )

    for i, (e5, w, g) in enumerate(zip(E4_BASE, W4_PAIR, g4), start=1):
        c52 = (A4 * w * R_52) * voigt_profile(E_ce, e5 + dE, sigma, g)
        c32 = (A4 * w * R_32) * voigt_profile(E_ce, e5 + DELTA_SO + dE, sigma, g)
        y += c52 + c32
        comps.append(
            {
                "spectrum": "Ce3d",
                "component": f"Ce4_pair{i}_3d5/2",
                "center_eV": float(e5 + dE),
                "sigma": sigma,
                "gamma": g,
                "fwhm_approx": voigt_fwhm_approx(sigma, g),
                "area": float(A4 * w * R_52),
            }
        )
        comps.append(
            {
                "spectrum": "Ce3d",
                "component": f"Ce4_pair{i}_3d3/2",
                "center_eV": float(e5 + DELTA_SO + dE),
                "sigma": sigma,
                "gamma": g,
                "fwhm_approx": voigt_fwhm_approx(sigma, g),
                "area": float(A4 * w * R_32),
            }
        )

    return y, comps


def gaussian_loglike_weighted(
    y: np.ndarray,
    ym: np.ndarray,
    log_a: float,
    log_b: float,
    confidence: np.ndarray,
) -> float:
    """
    Heteroscedastic Gaussian log-likelihood with confidence weighting.

    The unweighted variance is a + b * model intensity. The effective
    variance is inflated by 1 / C(E)^2 in low-confidence regions.
    """
    a = np.exp(log_a)
    b = np.exp(log_b)
    var = a + b * np.clip(ym, 0.0, None)
    var = np.clip(var, 1e-12, None)
    confidence = np.clip(confidence, 1e-6, 1.0)
    var_eff = var / (confidence * confidence)
    r = y - ym
    return float(-0.5 * np.sum((r * r) / var_eff + np.log(2.0 * np.pi * var_eff)))


def build_three_segment_tail_index(
    E_ce: np.ndarray,
    E_bg: np.ndarray,
    I_bg: np.ndarray,
) -> Tuple[np.ndarray, Dict[str, float | str]]:
    """
    Build the finalized Ce 3d scattering-tail accumulation index from a
    smart/Shirley background.

    Low/high flat windows are treated as full reset regions (T=0). Two
    partially flattened windows are treated as weak reset regions (T=0.25).
    The third segment starts from the Shirley background value at the right
    boundary of the weak-reset window to avoid an artificial discontinuity.
    """
    idx = np.argsort(E_bg)
    E_bg = E_bg[idx]
    I_bg = I_bg[idx]
    bg = np.interp(E_ce, E_bg, I_bg)

    platform_low = (870.78, 878.18)
    platform_mid = (891.98, 894.58)
    platform_weak = (910.0, 913.0)
    platform_high = (919.08, 924.58)
    weak_reset_t = 0.25

    def med(region: tuple[float, float]) -> float:
        lo, hi = min(region), max(region)
        m = (E_ce >= lo) & (E_ce <= hi)
        if np.sum(m) < 2:
            return float(np.median(bg))
        return float(np.median(bg[m]))

    def ramp(left: tuple[float, float], right: tuple[float, float], t0: float, t1: float) -> np.ndarray:
        b0 = med(left)
        b1 = med(right)
        denom = b1 - b0
        if abs(denom) < 1e-12:
            return np.full_like(E_ce, t0, dtype=float)
        r = np.clip((bg - b0) / denom, 0.0, 1.0)
        return t0 + (t1 - t0) * r

    def ramp_from_boundary(left_energy: float, right: tuple[float, float], t0: float, t1: float) -> np.ndarray:
        b0 = float(np.interp(left_energy, E_ce, bg))
        b1 = med(right)
        denom = b1 - b0
        if abs(denom) < 1e-12:
            return np.full_like(E_ce, t0, dtype=float)
        r = np.clip((bg - b0) / denom, 0.0, 1.0)
        return t0 + (t1 - t0) * r

    t = np.zeros_like(E_ce, dtype=float)
    low_start, low_end = platform_low
    mid_start, mid_end = platform_mid
    weak_start, weak_end = platform_weak
    high_start, high_end = platform_high

    seg1 = (E_ce >= low_end) & (E_ce <= mid_start)
    seg2 = (E_ce >= mid_end) & (E_ce <= weak_start)
    seg3 = (E_ce >= weak_end) & (E_ce <= high_start)
    t[seg1] = ramp(platform_low, platform_mid, 0.0, 1.0)[seg1]
    t[seg2] = ramp(platform_mid, platform_weak, weak_reset_t, 1.0)[seg2]
    t[seg3] = ramp_from_boundary(weak_end, platform_high, weak_reset_t, 1.0)[seg3]
    t[(E_ce >= low_start) & (E_ce <= low_end)] = 0.0
    t[(E_ce >= mid_start) & (E_ce <= mid_end)] = weak_reset_t
    t[(E_ce >= weak_start) & (E_ce <= weak_end)] = weak_reset_t
    t[(E_ce >= high_start) & (E_ce <= high_end)] = 0.0

    info: Dict[str, float | str] = {
        "B_low_platform": med(platform_low),
        "B_mid_platform": med(platform_mid),
        "B_weak_platform": med(platform_weak),
        "B_weak_right_boundary": float(np.interp(weak_end, E_ce, bg)),
        "B_high_platform": med(platform_high),
        "middle_weak_reset_tail_index": weak_reset_t,
        "weak_reset_tail_index": weak_reset_t,
        "third_segment_start": "right boundary of weak reset window",
    }
    return np.clip(t, 0.0, 1.0), info


def _normalize_01(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(values, dtype=float)
    return (values - lo) / (hi - lo)


def _contiguous_regions(mask: np.ndarray) -> list[tuple[int, int]]:
    idx = np.where(mask)[0]
    if idx.size == 0:
        return []
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends = np.r_[idx[breaks], idx[-1]]
    return [(int(s), int(e)) for s, e in zip(starts, ends)]


def detect_dynamic_confidence_starts(
    E_ce: np.ndarray,
    I_bg_sub: np.ndarray,
    *,
    smooth_window_eV: float = 1.6,
    envelope_threshold: float = 0.20,
    seed_threshold: float = 0.35,
    min_width_eV: float = 1.2,
    min_area_eV: float = 0.35,
    merge_gap_eV: float = 2.5,
    max_backtrack_eV: float = 4.0,
    min_positive_run_eV: float = 0.35,
) -> tuple[np.ndarray, Dict[str, float | str]]:
    """
    Detect dynamic confidence-reset positions from a background-subtracted Ce 3d spectrum.

    Major envelopes are selected with a strict hysteresis rule: a candidate
    envelope must cross the low normalized intensity threshold (0.20), contain a
    stronger seed (0.35), and satisfy minimum width/area requirements. For each
    retained envelope, the reset point is backtracked to the first point of the
    continuous positive-slope run before the envelope boundary.
    """
    order = np.argsort(E_ce)
    e = np.asarray(E_ce, dtype=float)[order]
    y = np.asarray(I_bg_sub, dtype=float)[order]
    if e.size < 12:
        return np.array([], dtype=float), {"dynamic_reset_count": 0.0}

    step = float(np.median(np.diff(e)))
    if not np.isfinite(step) or step <= 0:
        return np.array([], dtype=float), {"dynamic_reset_count": 0.0}
    win = int(round(smooth_window_eV / step))
    if win % 2 == 0:
        win += 1
    win = max(7, min(win, e.size - 1 if (e.size - 1) % 2 == 1 else e.size - 2))
    if win < 7:
        return np.array([], dtype=float), {"dynamic_reset_count": 0.0}

    y_s = savgol_filter(y, window_length=win, polyorder=3, mode="interp")
    y_norm = _normalize_01(y_s)
    d1 = np.gradient(y_s, e)

    active = (e >= 875.0) & (e <= 922.0)
    above = (y_norm >= envelope_threshold) & active
    regions = _contiguous_regions(above)

    merged_regions: list[tuple[int, int]] = []
    for start, end in regions:
        if not merged_regions:
            merged_regions.append((start, end))
            continue
        prev_start, prev_end = merged_regions[-1]
        if e[start] - e[prev_end] <= merge_gap_eV:
            merged_regions[-1] = (prev_start, end)
        else:
            merged_regions.append((start, end))

    positive_slope = d1 > 0
    min_positive_points = max(2, int(round(min_positive_run_eV / step)))
    backtrack_points = max(min_positive_points + 1, int(round(max_backtrack_eV / step)))

    def first_positive_run_before(region_start: int) -> int | None:
        search_start = max(0, region_start - backtrack_points)
        local_positive = positive_slope[search_start : region_start + 1]
        runs = _contiguous_regions(local_positive)
        candidates: list[tuple[int, int]] = []
        for start_local, end_local in runs:
            start = search_start + start_local
            end = search_start + end_local
            if end - start + 1 < min_positive_points:
                continue
            if y_norm[start] >= envelope_threshold:
                continue
            if end < region_start - min_positive_points:
                continue
            candidates.append((start, end))
        if not candidates:
            return None
        return candidates[-1][0]

    start_indices: list[int] = []
    retained_regions: list[tuple[int, int]] = []
    for start, end in merged_regions:
        width = float(e[end] - e[start])
        if width < min_width_eV:
            continue
        region_y = y_norm[start : end + 1]
        if float(np.max(region_y)) < seed_threshold:
            continue
        area = float(np.trapz(np.maximum(region_y - envelope_threshold, 0.0), e[start : end + 1]))
        if area < min_area_eV:
            continue
        onset = first_positive_run_before(start)
        if onset is None:
            continue
        start_indices.append(onset)
        retained_regions.append((start, end))

    start_energies = np.array([float(e[i]) for i in start_indices], dtype=float)
    info: Dict[str, float | str] = {
        "dynamic_reset_count": float(start_energies.size),
        "dynamic_reset_energies_eV": "; ".join(f"{v:.2f}" for v in start_energies),
        "dynamic_reset_rule": "v5 strict envelope seed, positive-slope backtracking",
        "envelope_threshold": float(envelope_threshold),
        "seed_threshold": float(seed_threshold),
        "min_width_eV": float(min_width_eV),
        "min_area_eV": float(min_area_eV),
        "merge_gap_eV": float(merge_gap_eV),
        "smooth_window_eV": float(smooth_window_eV),
        "smooth_window_points": float(win),
        "retained_envelope_regions_eV": "; ".join(f"{e[s]:.2f}-{e[t]:.2f}" for s, t in retained_regions),
    }
    return start_energies, info


def build_dynamic_confidence_tail_index(
    E_ce: np.ndarray,
    I_bg_sub: np.ndarray,
    E_bg: np.ndarray,
    I_bg: np.ndarray,
    *,
    weak_reset_t: float = 0.25,
) -> Tuple[np.ndarray, Dict[str, float | str]]:
    """
    Build the finalized dynamic confidence tail index.

    The first detected envelope onset is treated as the highest-confidence reset
    (T=0). All later detected onsets are weak resets (T=0.25). Between reset
    points, the index increases according to the local Shirley/background
    magnitude; if the background is locally flat, a linear energy ramp is used as
    a stable fallback.
    """
    order = np.argsort(E_ce)
    e = np.asarray(E_ce, dtype=float)[order]
    signal = np.asarray(I_bg_sub, dtype=float)[order]
    bg_order = np.argsort(E_bg)
    bg = np.interp(e, np.asarray(E_bg, dtype=float)[bg_order], np.asarray(I_bg, dtype=float)[bg_order])

    starts, info = detect_dynamic_confidence_starts(e, signal)
    t = np.zeros_like(e, dtype=float)
    if starts.size == 0:
        t_sorted = np.zeros_like(e, dtype=float)
    else:
        starts = np.sort(starts)
        reset_t = np.array([0.0] + [weak_reset_t] * (starts.size - 1), dtype=float)

        def ramp_between(mask: np.ndarray, left_e: float, right_e: float, left_t: float, right_t: float) -> None:
            if int(np.sum(mask)) == 0:
                return
            b0 = float(np.interp(left_e, e, bg))
            b1 = float(np.interp(right_e, e, bg))
            denom = b1 - b0
            if abs(denom) < 1e-12:
                denom_e = max(abs(right_e - left_e), 1e-12)
                r = np.clip((e[mask] - left_e) / denom_e, 0.0, 1.0)
            else:
                r = np.clip((bg[mask] - b0) / denom, 0.0, 1.0)
            t[mask] = left_t + (right_t - left_t) * r

        t[e <= starts[0]] = 0.0
        for i in range(starts.size - 1):
            mask = (e >= starts[i]) & (e <= starts[i + 1])
            ramp_between(mask, starts[i], starts[i + 1], reset_t[i], 1.0)
            reset_mask = np.isclose(e, starts[i + 1], atol=max(float(np.median(np.diff(e))) * 0.5, 1e-6))
            t[reset_mask] = reset_t[i + 1]
        mask = e >= starts[-1]
        ramp_between(mask, starts[-1], float(e[-1]), reset_t[-1], 1.0)
        t_sorted = np.clip(t, 0.0, 1.0)

    inverse = np.empty_like(order)
    inverse[order] = np.arange(order.size)
    out = t_sorted[inverse]
    info.update(
        {
            "weak_reset_tail_index": float(weak_reset_t),
            "first_reset_type": "highest confidence reset (T=0)",
            "later_reset_type": "weak reset (T=0.25)",
        }
    )
    return out, info
