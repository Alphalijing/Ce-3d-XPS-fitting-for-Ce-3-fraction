"""Short sampler, convergence-policy and output tests; no long data fit."""
from pathlib import Path
import inspect
import json
import sys
import tempfile
from threading import Event
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import emcee
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ce3d_xps_fitter_core as core
from ce3d_xps_model import Ce3dModel


def toy_proposal(ndim, one_branch=False):
    rng = np.random.default_rng(702)
    scales = np.geomspace(.1, 100, ndim)
    model = SimpleNamespace(ndim=ndim, scales=scales)
    train = rng.normal(size=(80, 40, ndim))
    train[..., 10] = np.abs(train[..., 10]) + .2 if one_branch else train[..., 10]
    pilot = train * scales
    best = np.zeros(ndim); best[10] = .7 * scales[10]
    alt = np.zeros(ndim); alt[10] = -.7 * scales[10]
    proposal, info = core.freeze_proposal(pilot, model, best, alt, .2 * scales)
    return model, proposal, info


class FrozenMixtureTests(unittest.TestCase):
    def test_complete_mixture_hastings_ratio_matches_scipy_11_and_15d(self):
        for ndim in (11, 15):
            for one_branch in (False, True):
                with self.subTest(ndim=ndim, fallback_branch=one_branch):
                    model, proposal, info = toy_proposal(ndim, one_branch)
                    rng = np.random.RandomState(734)
                    coords = rng.normal(size=(40, ndim)) * model.scales
                    new, correction = proposal(coords, rng)

                    def logq_physical(x):
                        terms = []
                        scale = np.diag(model.scales)
                        for w, mean, cov in zip(info['weights'], info['means_scaled'], info['covariances_scaled']):
                            terms.append(np.log(w) + multivariate_normal.logpdf(
                                x, mean=np.asarray(mean) * model.scales,
                                cov=scale @ np.asarray(cov) @ scale))
                        return logsumexp(terms, axis=0)

                    np.testing.assert_allclose(correction, logq_physical(coords) - logq_physical(new), rtol=2e-11, atol=2e-10)
                    self.assertTrue(info['frozen'] and info['full_mixture_Hastings_correction'])
                    self.assertTrue(info['untruncated_normal'])
                    self.assertEqual(sum(info['training_counts']), 40 * 40)
                    if one_branch:
                        self.assertEqual(info['training_counts'][1], 0)

    def test_reset_and_atomic_chunks_preserve_exact_production_trajectory(self):
        for ndim in (11, 15):
            with self.subTest(ndim=ndim), tempfile.TemporaryDirectory() as tmp:
                model, proposal, _ = toy_proposal(ndim)
                rng = np.random.RandomState(543)
                walkers = 40
                coords = rng.normal(size=(walkers, ndim))
                def logp(ts):
                    return -.5 * np.sum(np.asarray(ts)**2, axis=1)
                pilot = emcee.EnsembleSampler(walkers, ndim, logp, vectorize=True,
                    moves=[(emcee.moves.DEMove(), .8), (emcee.moves.DESnookerMove(), .2)])
                pilot.random_state = rng.get_state()
                state = pilot.run_mcmc(coords, 20)
                def production():
                    sampler = emcee.EnsembleSampler(walkers, ndim, logp, vectorize=True,
                        moves=[(emcee.moves.DEMove(), .5), (emcee.moves.DESnookerMove(), .2),
                               (emcee.moves.MHMove(proposal), .3)])
                    sampler.random_state = state.random_state
                    return sampler
                full = production()
                final_full = full.run_mcmc(emcee.State(state, copy=True), 160)
                chunks = production()
                cursor = emcee.State(state, copy=True)
                store = core.ChainStore(tmp)
                for length in (40, 60, 60):
                    cursor = chunks.run_mcmc(cursor, length)
                    store.append(chunks.get_chain(), chunks.get_log_prob(), cursor)
                    chunks.reset()
                    self.assertEqual(chunks.iteration, 0)
                actual, lp = store.read()
                np.testing.assert_array_equal(actual, full.get_chain())
                np.testing.assert_array_equal(lp, full.get_log_prob())
                np.testing.assert_array_equal(cursor.coords, final_full.coords)
                self.assertEqual(store.steps, 160)
                manifest = json.loads((Path(tmp) / 'chain' / 'manifest.json').read_text())
                self.assertEqual([(b['start'], b['end']) for b in manifest['blocks']], [(0, 40), (40, 100), (100, 160)])
                with np.load(Path(tmp) / 'chain' / 'last_state.npz') as saved:
                    np.testing.assert_array_equal(saved['coords'], final_full.coords)
                    np.testing.assert_array_equal(saved['log_prob'], final_full.log_prob)
                random_state = json.loads((Path(tmp) / 'chain' / 'random_state.json').read_text())
                np.testing.assert_array_equal(random_state['keys'], final_full.random_state[1])
                self.assertEqual(random_state['position'], final_full.random_state[2])
                self.assertFalse(list((Path(tmp) / 'chain').glob('*.tmp')))


class ConvergenceDiagnosticTests(unittest.TestCase):
    @staticmethod
    def iid(ndim):
        rng = np.random.default_rng(110 + ndim)
        chain = rng.normal(size=(10000, 40, ndim))
        chain[..., 2] = .2 + chain[..., 2] * .002
        lps = -.5 * np.sum(chain**2, axis=2)
        return chain, lps

    def test_short_and_unstable_are_rejected_stable_multivariate_iid_passes(self):
        for ndim in (11, 15):
            with self.subTest(ndim=ndim):
                chain, lps = self.iid(ndim)
                short = core.convergence_diagnostics(chain[:4001], lps[:4001], 4000, [])
                self.assertFalse(short['evaluated'])
                self.assertFalse(short['passed'])
                history = []
                for n in (6000, 8000, 10000):
                    diagnostic = core.convergence_diagnostics(chain[:n], lps[:n], 4000, history)
                    history.append(diagnostic)
                    if n < 10000:
                        self.assertFalse(diagnostic['passed'])
                self.assertTrue(diagnostic['passed'], diagnostic)
                self.assertEqual(len(diagnostic['all_tau']), ndim + 1)
                self.assertTrue(all(diagnostic['checks'].values()))
                unstable = chain.copy()
                unstable[4000:6000, :, 2] += .02
                unstable[8000:10000, :, 2] -= .02
                bad = core.convergence_diagnostics(unstable, lps, 4000, history[-2:])
                self.assertFalse(bad['passed'])
                self.assertFalse(bad['checks']['three_window_quantile_range_below_0p1_pp'])
                # CeIII can be stable while another parameter remains very
                # correlated. All dimensions and log-posterior are checked.
                rng = np.random.default_rng(19)
                increments = rng.normal(size=(10000, 40))
                chain[..., -1] = np.cumsum(increments, axis=0)
                bad_other = core.convergence_diagnostics(chain, lps, 4000, history[-2:])
                self.assertFalse(bad_other['passed'])
                self.assertFalse(bad_other['checks']['retained_length_at_least_100_max_tau'])


class FreshInitializationAndOutputTests(unittest.TestCase):
    @staticmethod
    def spectrum(satellites=False):
        energy = np.linspace(870.1, 924.9, 110)
        model = Ce3dModel(energy, 1000 + 500 * np.sin(energy), satellites=satellites)
        theta = model.initial()
        raw = model.predict(theta)[0] + 800
        return SimpleNamespace(energy=energy, raw=raw, background=np.full(len(energy), 800.),
            metadata={}, notes=[], source=ROOT / 'examples' / 'validation_sample_0p5.csv')

    def test_fresh_data_initialization_exercises_three_ls_and_four_map_starts(self):
        for satellites in (False, True):
            spectrum = self.spectrum(satellites)
            model = Ce3dModel(spectrum.energy, spectrum.raw - spectrum.background, satellites=satellites)
            def ls(fun, x, **kwargs):
                return SimpleNamespace(x=np.asarray(x), fun=fun(x), success=True)
            def mini(fun, x, **kwargs):
                self.assertTrue(np.isfinite(fun(x)))
                return SimpleNamespace(x=np.asarray(x), success=True, message='test stationary candidate')
            with patch.object(core, 'least_squares', side_effect=ls) as ls_mock, patch.object(core, 'minimize', side_effect=mini) as map_mock:
                rng, best, alt, candidates, widths, info = core.initialize(model, 42)
            self.assertEqual(ls_mock.call_count, 3)
            self.assertEqual(map_mock.call_count, 4)
            self.assertEqual(info['least_squares_starts'], 3)
            self.assertEqual(len(best), model.ndim)
            self.assertTrue(np.isfinite(model.log_prob(best)[0]))
            self.assertTrue(np.isfinite(model.log_prob(alt)[0]))
            walk = core.initialize_walkers(model, rng, best, alt, candidates, widths, 40)
            self.assertEqual(walk.shape, (40, model.ndim))
            self.assertTrue(np.isfinite(model.log_prob(walk)).all())
            self.assertEqual(len(ls_mock.call_args.args[1]), 13 if satellites else 9)

    def test_satellite_posterior_export_keeps_total_area_and_ce3_definition(self):
        spectrum = self.spectrum(True)
        model = Ce3dModel(spectrum.energy, spectrum.raw - spectrum.background, satellites=True)
        theta = model.initial()
        theta[11:13] = [.15, .2]
        chain = np.tile(theta, (30, 32, 1))
        chain[..., 2] += np.linspace(-.01, .01, 30)[:, None]
        store = SimpleNamespace(steps=30, read=lambda: (chain, np.zeros((30, 32))))
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            result = core._export(out, spectrum, model, store, 10, 'quick_complete',
                {'background': 'provided'}, {'evaluated': False, 'passed': False}, .1)
            peaks = pd.read_csv(out / 'peak_table.csv')
            samples = pd.read_csv(out / 'posterior_samples_ce3d.csv')
            self.assertEqual(samples.shape, (20 * 32, 15))
            self.assertEqual(list(samples.columns), list(model.param_names))
            self.assertEqual(len(peaks), 12)
            median = result['posterior_median']
            np.testing.assert_allclose(peaks['window_area'].sum(), median['A_ce_total'], rtol=2e-14)
            np.testing.assert_allclose(peaks['window_area'][:4].sum() / peaks['window_area'].sum(), median['x'], rtol=2e-14)
            expected = 100 * np.quantile(chain[10:, :, 2], [.025, .5, .975])
            np.testing.assert_allclose([result['Ce3_percent'][k] for k in ('q025', 'median', 'q975')], expected)
            self.assertIn('satellites included', result['Ce3_area_definition'])
            self.assertTrue(result['model']['satellites_enabled'])
            self.assertIn('sensitivity', result['model']['satellite_interpretation'])

    def test_pre_cancelled_run_never_claims_completion_or_starts_optimization(self):
        for mode in ('quick', 'convergence'):
            stop = Event(); stop.set()
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                with patch.object(core, 'load_spectrum', return_value=self.spectrum()), patch.object(core, '_save_background'), patch.object(core, 'least_squares', side_effect=AssertionError('optimization must not start')):
                    output = core.run_fit([ROOT / 'examples' / 'validation_sample_0p5.csv'], outdir=tmp, mode=mode, stop_event=stop)
                result = json.loads((output / 'details' / 'fit_result.json').read_text(encoding='utf-8'))
                self.assertEqual(result['status'], 'cancelled')
                self.assertFalse(result['posterior_available'])
                self.assertEqual(result['total_steps'], 0)
                self.assertFalse((output / 'initial_positions.npy').exists())


class _TinySampler:
    """Deterministic fake for testing the outer loop, with no MCMC work."""
    def __init__(self, walkers, ndim, logp, **kwargs):
        self.n = 0
        self.coords = None
        self.random_state = np.random.RandomState(1).get_state()
    def run_mcmc(self, state, steps, **kwargs):
        self.coords = state.coords if hasattr(state, 'coords') else np.asarray(state)
        self.n += steps
        return SimpleNamespace(coords=self.coords, log_prob=np.zeros(len(self.coords)), random_state=self.random_state)
    def get_chain(self):
        return np.broadcast_to(self.coords, (self.n,) + self.coords.shape)
    def get_log_prob(self):
        return np.zeros((self.n, len(self.coords)))
    def reset(self):
        self.n = 0


class _TinyStore:
    def __init__(self, folder):
        self.steps = 0
    def append(self, chain, lp, state):
        self.steps += len(chain)
        self.coords = state.coords
    def read(self):
        return np.broadcast_to(self.coords, (self.steps,) + self.coords.shape), np.zeros((self.steps, len(self.coords)))


class RunPolicyTests(unittest.TestCase):
    def test_api_defaults(self):
        signature = inspect.signature(core.run_fit)
        expected = dict(model='LA', satellites=False, mode='quick', walkers=80,
            steps=8000, burn=4000, satellite_gamma=.8, max_steps=None)
        for key, value in expected.items():
            self.assertEqual(signature.parameters[key].default, value)
        self.assertEqual(core.PILOT_STEPS, 2000)
        self.assertEqual(core.CHECK_EVERY, 2000)

    def run_policy(self, mode, pass_at=12000, cancel_at=None, max_steps=None):
        stop = Event()
        spectrum = FreshInitializationAndOutputTests.spectrum()
        model = Ce3dModel(spectrum.energy, spectrum.raw - spectrum.background)
        initial = model.initial()
        checkpoints = []
        result_holder = {}
        def diagnostic(chain, lps, burn, history):
            return dict(evaluated=True, passed=len(chain) >= pass_at, steps=len(chain))
        def on_progress(event):
            if event['event'] == 'checkpoint':
                checkpoints.append(event['steps'])
                if cancel_at is not None and event['steps'] >= cancel_at:
                    stop.set()
        def export(out, spectrum, target, store, burn, status, config, diagnostics, elapsed):
            result_holder.update(status=status, steps=store.steps, diagnostics=diagnostics)
            return {'posterior_available': False}
        init = (np.random.default_rng(1), initial, initial, [{'parameters': initial.tolist()}],
                model.scales * .001, {})
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(core, 'load_spectrum', return_value=spectrum), patch.object(core, '_save_background'), \
                 patch.object(core, 'initialize', return_value=init), \
                 patch.object(core, 'initialize_walkers', return_value=np.tile(initial, (22, 1))), \
                 patch.object(core, 'ChainStore', _TinyStore), patch.object(core.emcee, 'EnsembleSampler', _TinySampler), \
                 patch.object(core, 'freeze_proposal', return_value=(None, {'frozen': True})) as freeze, \
                 patch.object(core, 'convergence_diagnostics', side_effect=diagnostic), \
                 patch.object(core, '_export', side_effect=export):
                core.run_fit([ROOT / 'examples' / 'validation_sample_0p5.csv'], outdir=tmp, mode=mode,
                    walkers=22, steps=8000, stop_event=stop, max_steps=max_steps, progress=on_progress)
            self.assertEqual(freeze.call_count, 1)
        return result_holder, checkpoints

    def test_quick_ends_at_8000_and_convergence_continues_past_quick_budget(self):
        quick, checks = self.run_policy('quick')
        self.assertEqual(quick['status'], 'quick_complete')
        self.assertEqual(checks, [2000, 4000, 6000, 8000])
        convergence, checks = self.run_policy('convergence')
        self.assertEqual(convergence['status'], 'converged')
        self.assertEqual(checks[-1], 12000)
        self.assertGreater(convergence['steps'], 8000)

    def test_cancellation_precedes_completion_and_explicit_budget_is_honored(self):
        cancelled, _ = self.run_policy('convergence', pass_at=10000, cancel_at=10000)
        self.assertEqual(cancelled['status'], 'cancelled')
        cancelled_quick, _ = self.run_policy('quick', cancel_at=8000)
        self.assertEqual(cancelled_quick['status'], 'cancelled')
        budget, _ = self.run_policy('convergence', pass_at=12000, max_steps=8000)
        self.assertEqual(budget['status'], 'budget_exhausted')
        self.assertEqual(budget['steps'], 8000)


if __name__ == '__main__':
    unittest.main()
