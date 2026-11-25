"""Physics and frozen-target regressions, independent of GUI and sampler.

The reference helpers below are an explicit frozen transcription of the
validated mcmc_budget_targets.py target, not calls back into the new model.
The bundled real spectrum is used from examples; CE3D_VALIDATION_CURVES can select the original curves.csv.
"""
from pathlib import Path
import csv
import os
import sys
import unittest

import numpy as np
from scipy.fft import irfft, rfft, rfftfreq
from scipy.special import gammaln, voigt_profile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ce3d_xps_model import Ce3dModel


def reference_predict(ts, e, shape):
    centers = np.array([880.6, 899.2, 885.2, 903.8, 882.6, 901.2, 888.8, 907.4, 898.4, 917.])
    weights = np.array([2.9/8.3*.6, 2.9/8.3*.4, 5.4/8.3*.6, 5.4/8.3*.4,
                        3.9/10.2*.6, 3.9/10.2*.4, 2.5/10.2*.6, 2.5/10.2*.4,
                        3.8/10.2*.6, 3.8/10.2*.4])
    gi = np.array([4, 4, 5, 5, 6, 6, 7, 7, 8, 8])
    grid = (np.arange(8192) - 4096) * .025
    freq = rfftfreq(8192, .025)
    power = np.where(grid >= 0, .9, 2.)
    z = e[None, :, None] - centers[None, None, :] - ts[:, 0, None, None]
    p = voigt_profile(z, ts[:, 1, None, None], ts[:, gi][:, None, :])
    if shape == "LA":
        g = ts[:, 6, None]
        h = (1 + (grid[None, :] / g)**2)**(-power[None, :])
        h /= g * np.sqrt(np.pi) / 2 * (np.exp(gammaln(.4) - gammaln(.9)) + np.exp(gammaln(1.5) - gammaln(2.)))
        conv = irfft(rfft(h, axis=1) * np.exp(-2 * np.pi**2 * ts[:, 1, None]**2 * freq[None, :]**2), n=8192, axis=1)
        for i in range(len(ts)):
            for j in (4, 5):
                p[i, :, j] = np.interp(z[i, :, j], grid, conv[i])
    p /= np.trapz(p, e, axis=1)[:, None, :]
    amps = ts[:, 3, None] * weights[None, :] * np.where(np.arange(10)[None, :] < 4, ts[:, 2, None], 1 - ts[:, 2, None])
    return np.sum(p * amps[:, None, :], axis=2)


def reference_prior(t, net):
    # Exact original core hard limits and unnormalised prior terms.
    dE, sigma_ce, x, A_total = map(float, t[:4])
    gammas = [float(g) for g in t[4:9]]
    log_a_ce, log_b_ce = map(float, t[9:11])
    if not (-1.0 < dE < 1.0):
        return -np.inf
    if not (.05 < sigma_ce < 1.5):
        return -np.inf
    if not (0.0 < x < 1.0):
        return -np.inf
    if not (1.0 < A_total < max(float(np.max(net)) * 200.0, 1e3)):
        return -np.inf
    if any(not (.10 < g < 2.0) for g in gammas):
        return -np.inf
    if not (-20.0 < log_a_ce < 20.0 and -20.0 < log_b_ce < 20.0):
        return -np.inf
    lp = 0.0
    lp += -.5 * ((sigma_ce - .5) / .2)**2
    lp += -.5 * ((dE - 0.0) / .2)**2
    for g in gammas:
        lp += -.5 * ((g - .5) / .3)**2
    lp += -.5 * ((log_a_ce - np.log(max(np.var(net) * .05, 1e-6))) / 2.0)**2
    lp += -.5 * ((log_b_ce - np.log(1e-2)) / 2.0)**2
    return float(lp)


def reference_log_prob(ts, e, net, shape):
    prior = np.array([reference_prior(t, net) for t in ts])
    valid = np.isfinite(prior)
    ans = np.full(len(ts), -np.inf)
    if valid.any():
        p = reference_predict(ts[valid], e, shape)
        var = np.exp(ts[valid, 9, None]) + np.exp(ts[valid, 10, None]) * p
        ans[valid] = prior[valid] - .5 * np.sum(np.log(2 * np.pi * var) + (net - p)**2 / var, axis=1)
    return ans


def real_spectrum():
    paths = []
    if os.environ.get("CE3D_VALIDATION_CURVES"):
        paths.append(Path(os.environ["CE3D_VALIDATION_CURVES"]))
    paths += list((ROOT / "examples").glob("*.csv"))
    for p in paths:
        if not p.is_file():
            continue
        with p.open(encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.DictReader(fh))
        if not rows or "energy_eV" not in rows[0]:
            continue
        if "raw" in rows[0] and "LA_fixed_background" in rows[0]:
            return np.array([float(r["energy_eV"]) for r in rows]), np.array([
                float(r["raw"]) - float(r["LA_fixed_background"]) for r in rows
            ])
        if "net" in rows[0]:
            return np.array([float(r["energy_eV"]) for r in rows]), np.array([float(r["net"]) for r in rows])
        if "raw" in rows[0] and "background" in rows[0]:
            return np.array([float(r["energy_eV"]) for r in rows]), np.array([
                float(r["raw"]) - float(r["background"]) for r in rows
            ])
    return None


def parameter_vectors(model):
    rng = np.random.default_rng(23419)
    ts = np.tile(model.initial(), (16, 1))
    ts[:, 0] = rng.uniform(-.98, .98, 16)
    ts[:, 1] = rng.uniform(.06, 1.49, 16)
    ts[:, 2] = rng.uniform(.01, .99, 16)
    ts[:, 3] *= rng.uniform(.3, 1.2, 16)
    ts[:, 3] = np.clip(ts[:, 3], 2, model.bounds[1][3] * .95)
    ts[:, 4:9] = rng.uniform(.101, 1.999, (16, 5))
    ts[:, 9] = rng.uniform(-19.0, 19.0, 16)
    ts[:, 10] = rng.uniform(-19.0, 19.0, 16)
    return ts


class FrozenTargetTests(unittest.TestCase):
    def compare_target(self, e, net):
        for shape in ("LA", "Voigt"):
            with self.subTest(shape=shape):
                model = Ce3dModel(e, net, lineshape=shape)
                ts = parameter_vectors(model)
                # Bitwise equality is expected: the numerical operations in
                # the original target have deliberately been preserved.
                np.testing.assert_array_equal(model.predict(ts), reference_predict(ts, e, shape))
                np.testing.assert_array_equal(model.log_prob(ts), reference_log_prob(ts, e, net, shape))
                np.testing.assert_array_equal([model.log_prior(t) for t in ts], [reference_prior(t, net) for t in ts])

    def test_sixteen_vectors_frozen_la_and_voigt(self):
        e = np.linspace(870, 925, 550)
        net = 1e4 * (1 + np.cos(e * .2)) + 500
        self.compare_target(e, net)

    def test_real_data_sixteen_vectors_frozen_la_and_voigt(self):
        data = real_spectrum()
        if data is None:
            self.skipTest("original validation spectrum is unavailable; set CE3D_VALIDATION_CURVES")
        self.compare_target(*data)

    def test_unchanged_open_prior_boundaries(self):
        e = np.linspace(870, 925, 550)
        net = 1e4 * (1 + np.cos(e * .2))
        model = Ce3dModel(e, net)
        for j in range(11):
            for edge in model.bounds:
                t = model.initial()
                t[j] = edge[j]
                with self.subTest(parameter=j, edge=edge[j]):
                    self.assertEqual(model.log_prior(t), -np.inf)
                    self.assertEqual(model.log_prior(t), reference_prior(t, net))
        t = model.initial()
        t[0] = -.999999999
        t[6] = .100000001
        self.assertTrue(np.isfinite(model.log_prob(t)[0]))


class SatelliteAndAreaTests(unittest.TestCase):
    def setUp(self):
        self.e = np.linspace(870.08, 924.98, 550)
        self.net = 2e4 * (1 + .7 * np.sin(self.e * .15))

    def test_observed_window_areas_and_ce3_fraction(self):
        for shape in ("LA", "Voigt"):
            for sat in (False, True):
                with self.subTest(shape=shape, satellites=sat):
                    model = Ce3dModel(self.e, self.net, lineshape=shape, satellites=sat)
                    t = model.initial()
                    t[0], t[1], t[2] = -.8, 1.2, .37
                    if sat:
                        t[11:13] = [.18, .27]
                    c = model.components(t)
                    areas = np.trapz(c, self.e, axis=0)
                    self.assertEqual(c.shape, (len(self.e), 12 if sat else 10))
                    np.testing.assert_allclose(c.sum(axis=1), model.predict(t)[0], rtol=1e-15)
                    np.testing.assert_allclose(areas.sum(), t[3], rtol=1e-14)
                    np.testing.assert_allclose(areas[:4].sum() / areas.sum(), t[2], rtol=1e-14)
                    np.testing.assert_array_equal(model.ce3_fraction(t), [t[2]])
                    np.testing.assert_allclose([r["window_area"] for r in model.peak_table(t)], areas, rtol=1e-15)
                    if sat:
                        np.testing.assert_allclose(areas[10:], t[3] * (1 - t[2]) * t[11:13], rtol=1e-14)
                        np.testing.assert_allclose(areas[4:10].sum(), t[3] * (1 - t[2]) * (1 - t[11:13].sum()), rtol=1e-14)
                    # Every template doublet retains the fixed 3:2 branches.
                    np.testing.assert_allclose(areas[:10:2] / areas[1:10:2], 1.5, rtol=1e-14)

    def test_satellites_zero_recovers_off_exactly(self):
        for shape in ("LA", "Voigt"):
            off = Ce3dModel(self.e, self.net, lineshape=shape)
            on = Ce3dModel(self.e, self.net, lineshape=shape, satellites=True)
            t = off.initial()
            s = np.r_[t, 0, 0, 897.5, 914.75]
            np.testing.assert_array_equal(off.predict(t), on.predict(s))
            np.testing.assert_array_equal(off.log_prob(t), on.log_prob(s))
            self.assertEqual(off.log_prior(t), on.log_prior(s))

    def test_satellite_centers_are_absolute_and_shape_is_voigt(self):
        model = Ce3dModel(self.e, self.net, satellites=True, satellite_gamma=.9)
        t = model.initial()
        s = t.copy()
        t[0], s[0] = -.7, .8
        np.testing.assert_array_equal(model.components(t)[:, 10:], model.components(s)[:, 10:])
        a, b = model.peak_table(t), model.peak_table(s)
        np.testing.assert_array_equal([r["center_eV"] for r in a[10:]], [r["center_eV"] for r in b[10:]])
        for r in a[10:]:
            self.assertEqual(r["lineshape"], "Voigt")
            self.assertEqual(r["width_parameter_gamma_eV"], .9)
            self.assertEqual(r["Lorentzian_HWHM_eV"], .9)
            self.assertIsNone(r["LA_kernel_gamma_eV"])
            self.assertIn("sensitivity", r["interpretation"])
        for j in range(10):
            self.assertAlmostEqual(b[j]["center_eV"] - a[j]["center_eV"], 1.5)
        for j in range(0, 10, 2):
            self.assertAlmostEqual(a[j + 1]["center_eV"] - a[j]["center_eV"], 18.6)

    def test_gamma_columns_distinguish_la_scale_from_lorentzian_hwhm(self):
        for shape in ("LA", "Voigt"):
            model = Ce3dModel(self.e, self.net, lineshape=shape, satellites=True)
            t = model.initial()
            table = model.peak_table(t)
            for j, row in enumerate(table):
                is_la = shape == "LA" and j in (4, 5)
                with self.subTest(shape=shape, component=j):
                    gamma = row["width_parameter_gamma_eV"]
                    self.assertNotIn("gamma_HWHM_eV", row)
                    if is_la:
                        self.assertEqual(gamma, t[6])
                        self.assertEqual(row["LA_kernel_gamma_eV"], gamma)
                        self.assertIsNone(row["Lorentzian_HWHM_eV"])
                        self.assertIn("not the kernel or convolved-peak HWHM", row["gamma_definition"])
                        self.assertIn("g*sqrt(2**(1/power)-1)", row["gamma_definition"])
                    else:
                        self.assertEqual(row["Lorentzian_HWHM_eV"], gamma)
                        self.assertIsNone(row["LA_kernel_gamma_eV"])
                        self.assertIn("Lorentzian kernel HWHM", row["gamma_definition"])
                    if j >= 10:
                        self.assertEqual(gamma, .8)
            metadata = model.metadata()
            self.assertIn("LA", metadata["gamma_definition"])
            self.assertIn("Voigt", metadata["gamma_definition"])
            self.assertEqual(metadata["gamma_definition"]["g4_1"], table[4]["gamma_definition"])
            self.assertEqual(metadata["gamma_definition"]["satellite_gamma"], table[10]["gamma_definition"])

    def test_satellite_fraction_simplex_is_uniform_and_nonnegative(self):
        model = Ce3dModel(self.e, self.net, satellites=True)
        t = model.initial()
        for fractions in ([0, 0], [.7, .2], [.05, .9]):
            s = t.copy()
            s[11:13] = fractions
            self.assertEqual(model.log_prior(s), model.log_prior(t))
        for fractions in ([-.001, .1], [.1, -.001], [.6, .4], [.7, .5], [1.1, 0]):
            s = t.copy()
            s[11:13] = fractions
            with self.subTest(fractions=fractions):
                self.assertEqual(model.log_prior(s), -np.inf)
                self.assertEqual(model.log_prob(s)[0], -np.inf)

    def test_center_bounds_and_finite_initialization(self):
        for net in (self.net, np.zeros_like(self.e), -np.ones_like(self.e)):
            for sat in (False, True):
                model = Ce3dModel(self.e, net, satellites=sat)
                self.assertEqual(model.ndim, 15 if sat else 11)
                self.assertEqual(len(model.param_names), model.ndim)
                self.assertEqual(model.scales.shape, (model.ndim,))
                self.assertTrue(np.all(np.isfinite(model.scales)))
                self.assertTrue(np.all(model.scales > 0))
                lo, hi = model.bounds
                self.assertTrue(np.all(np.isfinite(lo)) and np.all(np.isfinite(hi)))
                self.assertTrue(np.all(hi > lo))
                self.assertTrue(np.isfinite(model.log_prob(model.initial())[0]))
        model = Ce3dModel(self.e, self.net, satellites=True)
        for j, edges in ((13, [896.999, 898.001]), (14, [914.499, 915.001])):
            for value in edges:
                s = model.initial()
                s[j] = value
                self.assertEqual(model.log_prior(s), -np.inf)

    def test_only_ce4_v_u_change_with_lineshape(self):
        a = Ce3dModel(self.e, self.net)
        b = Ce3dModel(self.e, self.net, lineshape="Voigt")
        t = a.initial()
        ca, cb = a.components(t), b.components(t)
        np.testing.assert_array_equal(ca[:, [0, 1, 2, 3, 6, 7, 8, 9]], cb[:, [0, 1, 2, 3, 6, 7, 8, 9]])
        self.assertGreater(np.max(np.abs(ca[:, 4:6] - cb[:, 4:6])), 1)

    def test_input_errors(self):
        for shape in ("Gaussian", ""):
            with self.assertRaises(ValueError):
                Ce3dModel(self.e, self.net, lineshape=shape)
        for gamma in (0, -1, np.nan, np.inf):
            with self.assertRaises(ValueError):
                Ce3dModel(self.e, self.net, satellite_gamma=gamma)
        for e in (self.e[::-1], np.full_like(self.e, 900), self.e * np.nan):
            with self.assertRaises(ValueError):
                Ce3dModel(e, self.net)


if __name__ == "__main__":
    unittest.main()
