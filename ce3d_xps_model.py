"""Ce 3d observed-window area model and the validated Bayesian target.

LA is the default: alpha=0.9 and beta=2 only for Ce(IV) v/u. All other
main peaks are Voigt. The first eleven coordinates and their prior exactly
match the validated target. Optional satellites are a sensitivity model;
their intensities have no externally calibrated reference ratio.
"""
from __future__ import annotations

import numpy as np
from scipy.fft import irfft, rfft, rfftfreq
from scipy.special import gammaln, voigt_profile


CE_COLS = [
    "dE", "sigma_ce", "x", "A_ce_total", "g3_1", "g3_2", "g4_1",
    "g4_2", "g4_3", "log_a_ce", "log_b_ce",
]
SATELLITE_COLS = [
    "satellite_1_fraction_Ce4", "satellite_2_fraction_Ce4",
    "satellite_1_center_eV", "satellite_2_center_eV",
]
MAIN_CENTERS = np.array([
    880.6, 899.2, 885.2, 903.8, 882.6, 901.2, 888.8, 907.4, 898.4, 917.0,
])
MAIN_WEIGHTS = np.array([
    2.9/8.3*.6, 2.9/8.3*.4, 5.4/8.3*.6, 5.4/8.3*.4,
    3.9/10.2*.6, 3.9/10.2*.4, 2.5/10.2*.6, 2.5/10.2*.4,
    3.8/10.2*.6, 3.8/10.2*.4,
])
GAMMA_INDICES = np.array([4, 4, 5, 5, 6, 6, 7, 7, 8, 8])
MAIN_COMPONENT_NAMES = (
    "Ce3+ v0", "Ce3+ u0", "Ce3+ v'", "Ce3+ u'", "Ce4+ v", "Ce4+ u",
    "Ce4+ v''", "Ce4+ u''", "Ce4+ v'''", "Ce4+ u'''",
)
LA_ALPHA = 0.9
LA_BETA = 2.0
LA_GRID_N = 8192
LA_GRID_DX = 0.025
VOIGT_GAMMA_DEFINITION = (
    "Lorentzian kernel HWHM before Gaussian convolution; not the total Voigt peak HWHM"
)
LA_GAMMA_DEFINITION = (
    "Powered-Lorentz LA kernel scale g in (1+(z/g)^2)^(-power); "
    "side-specific kernel HWHM is g*sqrt(2**(1/power)-1), with power=0.9 "
    "at z>=0 and 2 at z<0; g is not the kernel or convolved-peak HWHM"
)


class Ce3dModel:
    """Fit a fixed-background net spectrum on increasing binding energy.

    ``A_ce_total`` is the integrated total signal over the observed energy
    window, including satellites when enabled. ``x`` is consequently the
    Ce(III) fraction of that same total window area in both model variants.
    Satellite coordinates are absolute energies, unaffected by ``dE``.
    """

    def __init__(self, energy, net, lineshape="LA", satellites=False,
                 satellite_gamma=0.8):
        self.energy = np.asarray(energy, dtype=float).copy()
        self.net = np.asarray(net, dtype=float).copy()
        if (self.energy.ndim != 1 or self.net.shape != self.energy.shape
                or self.energy.size < 2):
            raise ValueError("energy and net must be matching one-dimensional arrays with at least two points")
        if not (np.all(np.isfinite(self.energy)) and np.all(np.isfinite(self.net))):
            raise ValueError("energy and net must contain only finite values")
        if not np.all(np.diff(self.energy) > 0):
            raise ValueError("energy must be strictly increasing; sort the spectrum before fitting")
        shapes = {"la": "LA", "voigt": "Voigt"}
        if str(lineshape).lower() not in shapes:
            raise ValueError("lineshape must be 'LA' or 'Voigt'")
        self.lineshape = shapes[str(lineshape).lower()]
        self.satellites = bool(satellites)
        self.satellite_gamma = float(satellite_gamma)
        if not (np.isfinite(self.satellite_gamma) and self.satellite_gamma > 0):
            raise ValueError("satellite_gamma must be finite and positive (Lorentzian HWHM, eV)")
        self.param_names = tuple(CE_COLS + (SATELLITE_COLS if self.satellites else []))
        self.ndim = len(self.param_names)
        self.component_names = MAIN_COMPONENT_NAMES + (
            ("Ce4+ satellite 1", "Ce4+ satellite 2") if self.satellites else ()
        )
        self.n_components = len(self.component_names)
        self._amp_max = max(float(np.max(self.net)) * 200.0, 1e3)
        self._log_a_mean = np.log(max(np.var(self.net) * 0.05, 1e-6))
        lo = np.r_[[-1.0, .05, 0.0, 1.0], [.1] * 5, -20.0, -20.0]
        hi = np.r_[[1.0, 1.5, 1.0, self._amp_max], [2.0] * 5, 20.0, 20.0]
        if self.satellites:
            lo = np.r_[lo, 0.0, 0.0, 897.0, 914.5]
            hi = np.r_[hi, 1.0, 1.0, 898.0, 915.0]
        self.bounds = (lo, hi)
        initial_area = max(float(np.trapz(np.clip(self.net, 0, None), self.energy)), 100.0)
        # Numerical coordinates only: these scales never enter the posterior.
        area_scale = np.clip(initial_area, 100.0, self._amp_max)
        self.scales = np.array([.2, .2, .1, area_scale, .3, .3, .3, .3, .3, 2., 2.])
        if self.satellites:
            self.scales = np.r_[self.scales, .03, .03, .2, .1]
        self._grid = (np.arange(LA_GRID_N) - LA_GRID_N // 2) * LA_GRID_DX
        self._freq = rfftfreq(LA_GRID_N, LA_GRID_DX)
        self._power = np.where(self._grid >= 0, LA_ALPHA, LA_BETA)

    def initial(self):
        """Return a finite starting point derived directly from this spectrum."""
        area = max(float(np.trapz(np.clip(self.net, 0, None), self.energy)), 100.0)
        theta = np.array([
            0.0, .55, .2, area, .45, .55, .45, .55, .65,
            self._log_a_mean, np.log(1e-2),
        ])
        if self.satellites:
            theta = np.r_[theta, .01, .01, 897.5, 914.75]
        lo, hi = self.bounds
        margin = np.maximum((hi - lo) * 1e-9, 1e-9)
        # Ordinary spectra reproduce the original initializer. Clipping only
        # handles unusual dynamic ranges that place the raw seed out of prior.
        return np.clip(theta, lo + margin, hi - margin)

    def _batch(self, theta):
        ts = np.asarray(theta, dtype=float)
        if ts.ndim == 1:
            ts = ts[None, :]
        if ts.ndim != 2 or ts.shape[1] != self.ndim:
            raise ValueError(f"expected parameter array with final dimension {self.ndim}")
        return ts

    def _unit_profiles(self, ts):
        e = self.energy
        z = e[None, :, None] - MAIN_CENTERS[None, None, :] - ts[:, 0, None, None]
        p = voigt_profile(z, ts[:, 1, None, None], ts[:, GAMMA_INDICES][:, None, :])
        if self.lineshape == "LA":
            g = ts[:, 6, None]
            h = (1 + (self._grid[None, :] / g) ** 2) ** (-self._power[None, :])
            h /= g * np.sqrt(np.pi) / 2 * (
                np.exp(gammaln(.4) - gammaln(.9))
                + np.exp(gammaln(1.5) - gammaln(2.))
            )
            conv = irfft(rfft(h, axis=1) * np.exp(
                -2 * np.pi**2 * ts[:, 1, None]**2 * self._freq[None, :]**2
            ), n=LA_GRID_N, axis=1)
            for i in range(len(ts)):
                for j in (4, 5):
                    p[i, :, j] = np.interp(z[i, :, j], self._grid, conv[i])
        p /= np.trapz(p, e, axis=1)[:, None, :]
        if self.satellites:
            sat_z = e[None, :, None] - ts[:, None, 13:15]
            sat = voigt_profile(sat_z, ts[:, 1, None, None], self.satellite_gamma)
            sat /= np.trapz(sat, e, axis=1)[:, None, :]
            p = np.concatenate((p, sat), axis=2)
        return p

    def _amplitudes(self, ts):
        amps = ts[:, 3, None] * MAIN_WEIGHTS[None, :] * np.where(
            np.arange(10)[None, :] < 4, ts[:, 2, None], 1 - ts[:, 2, None]
        )
        if self.satellites:
            f = ts[:, 11:13]
            amps[:, 4:] *= 1 - f.sum(axis=1)[:, None]
            sat_amps = ts[:, 3, None] * (1 - ts[:, 2, None]) * f
            amps = np.concatenate((amps, sat_amps), axis=1)
        return amps

    def predict(self, theta_batch):
        """Return peak signal, shape (number of parameter vectors, energies)."""
        ts = self._batch(theta_batch)
        return np.sum(self._unit_profiles(ts) * self._amplitudes(ts)[:, None, :], axis=2)

    def log_prior(self, theta):
        """Original core prior plus optional uniform satellite simplex/bounds.

        The satellite fractions are jointly uniform on f1>=0, f2>=0,
        f1+f2<1. No ratio penalty or reference intensity is imposed. As in
        the original core, parameter-independent prior constants are omitted.
        """
        t = np.asarray(theta, dtype=float)
        if t.shape != (self.ndim,) or not np.all(np.isfinite(t)):
            return -np.inf
        dE, sigma_ce, x, A_total = map(float, t[:4])
        gammas = [float(g) for g in t[4:9]]
        log_a_ce, log_b_ce = map(float, t[9:11])
        if not (-1.0 < dE < 1.0):
            return -np.inf
        if not (0.05 < sigma_ce < 1.5):
            return -np.inf
        if not (0.0 < x < 1.0):
            return -np.inf
        if not (1.0 < A_total < self._amp_max):
            return -np.inf
        if any(not (0.10 < g < 2.0) for g in gammas):
            return -np.inf
        if not (-20.0 < log_a_ce < 20.0 and -20.0 < log_b_ce < 20.0):
            return -np.inf
        if self.satellites:
            f1, f2, c1, c2 = t[11:15]
            if not (f1 >= 0.0 and f2 >= 0.0 and f1 + f2 < 1.0
                    and 897.0 <= c1 <= 898.0 and 914.5 <= c2 <= 915.0):
                return -np.inf
        lp = 0.0
        lp += -0.5 * ((sigma_ce - 0.5) / 0.2) ** 2
        lp += -0.5 * ((dE - 0.0) / 0.2) ** 2
        for g in gammas:
            lp += -0.5 * ((g - 0.5) / 0.3) ** 2
        lp += -0.5 * ((log_a_ce - np.log(max(np.var(self.net) * 0.05, 1e-6))) / 2.0) ** 2
        lp += -0.5 * ((log_b_ce - np.log(1e-2)) / 2.0) ** 2
        return float(lp)

    def log_prob(self, theta_batch):
        """Original prior and Gaussian likelihood with var=a+b*peak_signal."""
        ts = self._batch(theta_batch)
        prior = np.array([self.log_prior(t) for t in ts])
        valid = np.isfinite(prior)
        ans = np.full(len(ts), -np.inf)
        if valid.any():
            p = self.predict(ts[valid])
            var = np.exp(ts[valid, 9, None]) + np.exp(ts[valid, 10, None]) * p
            ans[valid] = prior[valid] - .5 * np.sum(
                np.log(2 * np.pi * var) + (self.net - p) ** 2 / var, axis=1
            )
        return ans

    def components(self, theta):
        """Actual area-scaled component curves, shape (energies, components)."""
        t = np.asarray(theta, dtype=float)
        if t.shape != (self.ndim,):
            raise ValueError(f"components expects one {self.ndim}-parameter vector")
        ts = t[None, :]
        return (self._unit_profiles(ts) * self._amplitudes(ts)[:, None, :])[0]

    def ce3_fraction(self, theta_batch):
        """Ce(III) fraction of total modeled area in the observed window."""
        return self._batch(theta_batch)[:, 2].copy()

    def peak_table(self, theta):
        """Return measured-window component areas in intensity-times-eV units."""
        t = np.asarray(theta, dtype=float)
        curves = self.components(t)
        areas = np.trapz(curves, self.energy, axis=0)
        centers = MAIN_CENTERS + t[0]
        gammas = t[GAMMA_INDICES]
        if self.satellites:
            centers = np.r_[centers, t[13:15]]
            gammas = np.r_[gammas, self.satellite_gamma, self.satellite_gamma]
        rows = []
        for j, name in enumerate(self.component_names):
            satellite = j >= 10
            is_la = self.lineshape == "LA" and j in (4, 5)
            rows.append({
                "component": name,
                "oxidation_state": "Ce(III)" if j < 4 else "Ce(IV)",
                "center_eV": float(centers[j]),
                "sigma_eV": float(t[1]),
                "width_parameter_gamma_eV": float(gammas[j]),
                "Lorentzian_HWHM_eV": None if is_la else float(gammas[j]),
                "LA_kernel_gamma_eV": float(gammas[j]) if is_la else None,
                "gamma_definition": LA_GAMMA_DEFINITION if is_la else VOIGT_GAMMA_DEFINITION,
                "lineshape": "LA" if is_la else "Voigt",
                "area": float(areas[j]),
                "window_area": float(areas[j]),
                "fraction_of_total": float(areas[j] / t[3]),
                "satellite": satellite,
                "interpretation": "sensitivity model; no independently calibrated intensity ratio" if satellite else "fixed main-peak template",
            })
        return rows

    def metadata(self):
        """Machine-readable model assumptions for fit reports."""
        return {
            "lineshape": self.lineshape,
            "parameter_count": self.ndim,
            "parameter_names": list(self.param_names),
            "component_normalization": "each component integrates to one over the observed energy window before area scaling",
            "energy_window_eV": [float(self.energy[0]), float(self.energy[-1])],
            "doublet_split_eV": 18.6,
            "doublet_branch_ratio": "3:2",
            "main_centers_eV": MAIN_CENTERS.tolist(),
            "main_branch_weights": MAIN_WEIGHTS.tolist(),
            "shared_sigma": True,
            "gamma_definition": {
                "Voigt": VOIGT_GAMMA_DEFINITION,
                "LA": LA_GAMMA_DEFINITION,
                "g4_1": LA_GAMMA_DEFINITION if self.lineshape == "LA" else VOIGT_GAMMA_DEFINITION,
                "satellite_gamma": VOIGT_GAMMA_DEFINITION if self.satellites else None,
            },
            "LA_alpha": LA_ALPHA,
            "LA_beta": LA_BETA,
            "LA_components": [4, 5] if self.lineshape == "LA" else [],
            "LA_grid_points": LA_GRID_N,
            "LA_grid_spacing_eV": LA_GRID_DX,
            "variance": "exp(log_a_ce)+exp(log_b_ce)*peak_model",
            "prior": "original eleven-parameter core prior; unchanged hard limits and Gaussian terms",
            "Ce3_fraction_definition": "Ce(III) observed-window area / total Ce observed-window area including satellites",
            "satellites_enabled": self.satellites,
            "satellite_interpretation": "sensitivity model; no independently calibrated reference intensity ratio" if self.satellites else None,
            "satellite_fraction_prior": "uniform simplex: f1>=0, f2>=0, f1+f2<1" if self.satellites else None,
            "satellite_center_ranges_eV": [[897.0, 898.0], [914.5, 915.0]] if self.satellites else [],
            "satellite_centers_absolute": self.satellites,
            "satellite_gamma_HWHM_eV": self.satellite_gamma if self.satellites else None,
        }
