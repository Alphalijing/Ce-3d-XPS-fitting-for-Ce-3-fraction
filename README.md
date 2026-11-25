# Ce 3d XPS Fitter 2.0.1

Automated, physics-constrained Bayesian quantification of the Ce³⁺ fraction from Ce 3d XPS spectra. The default workflow uses **Gaussian-convolved LA-type peaks with fixed α = 0.9 and β = 2**, a **fixed Shirley background**, and **MCMC posterior sampling**.

## Start fitting

**Windows:** download `Ce3d_XPS_Fitter_v2.0.1_Windows_EN.zip` from this repository's [Releases](https://github.com/Alphalijing/Ce-3d-XPS-fitting-for-Ce-3-fraction/releases), extract it completely, and open `Ce3d_XPS_Fitter_v2.0.1_EN.exe`. Keep its adjacent `_internal` folder.

**Python source:** Python 3.11 is recommended.

```console
python -m pip install -r requirements.txt
python ce3d_xps_fitter_gui.py
```

The interface provides three independent choices:

| Choice | Default | Alternative |
|---|---|---|
| Peak model | LA-type | Voigt |
| Satellite contributions | Exclude | Include |
| Sampling mode | Quick fit: 8000 steps | Convergence validation |

The default is 80 walkers, 8000 total steps including a 2000-step pilot, with the first 4000 steps discarded. Quick fit performs MCMC. Progress is reported at the start and every 2000 steps. Convergence validation continues until the implemented diagnostic criteria pass or the user stops it.

Supported inputs include instrument `.xls/.xlsx`, two-column text/CSV, and regular XPS scans in `.vms/.vamas`. Three-column text/CSV can supply an explicitly fixed background. Multiple VAMAS blocks/channels require selection; kinetic-energy axes require an explicit energy offset. See the [interface guide](GUI_GUIDE.md) and [application reference](docs/APPLICATION_REFERENCE.md).

## Results

`Results.html` presents **Ce³⁺ % ± one posterior standard deviation** first, followed by component fits, concise experimental conditions, peak parameters and posterior distributions. The median and 95% equal-tail credible interval are reported separately. Full records include `details/experiment_metadata.json`, which preserves available acquisition conditions, including acquisition time where supplied, selected scan/channel and source hash. Missing conditions remain unspecified.

## Publication data and reproduction

The repository contains raw inputs, frozen backgrounds, fitted curves, parameter summaries and reproduction instructions for the revised manuscript.

| Dataset | Contents |
|---|---|
| [Six La-doped specimens](data/experimental/la_series) | Original Excel exports; composition/catalysis data; LA and Voigt MCMC summaries; DS point estimates; all Ce³⁺ posterior draws used for the reported intervals |
| [19 additional spectra](data/experimental/additional_19) | Raw Excel exports, saved background constructions, LA/Voigt fits and high-confidence-region residuals |
| [48 synthetic cases](data/synthetic) | 24 background-free and 24 cumulative-tail cases; exact truth and noise seeds; matched constrained-MCMC and independent-parameter-LA fits |
| [Further validation](data/validation) | Background comparisons, empirical Tougaard-form construction, satellites, joint backgrounds, sampling budgets and prior sensitivity |
| [Input examples](examples) | Excel, CSV and an attributed public multi-block VAMAS example |

Main-text Figure 4 uses LA and Voigt posterior medians and 95% credible intervals; DS retains the deterministic comparison. The 19-spectrum extension retains its archived deterministic fits. These estimators are labelled in the data and [reproduction guide](docs/REPRODUCIBILITY.md).

From the repository root:

```console
python analysis/verify_release.py
python analysis/summarize_results.py
python analysis/plot_figure3.py
python analysis/plot_figure4.py
python analysis/plot_figure5.py
python analysis/plot_posterior_s2.py
```

These commands use saved data and **do not start MCMC**. Figures and recalculated statistics are written under `results/`, leaving the published data intact. The 19-spectrum extension is supplied as tables and data without a separate figure.

For the exact benchmark protocol, data mapping, synthetic generator and optional refitting commands, see [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md). For scientific scope and mathematical definitions, see [MODEL_AND_SCOPE.md](docs/MODEL_AND_SCOPE.md).

## Tests and build

```console
python -m pip install -r requirements-dev.txt
python -m pytest tests -q
python build_windows_exe.py
```

The production scientific source is unchanged from the validated English v2.0.1 release. This publication bundle adds datasets and portable analysis scripts. See [VALIDATION.md](VALIDATION.md) for inherited checks and [publication validation](docs/PUBLICATION_VALIDATION.md) for the checks performed on this bundle.

## Interpretation and citation

Posterior intervals are conditional on the selected background, profiles, constraints, priors and noise model. They do not represent all experimental or model uncertainty. The LA implementation does not claim identity with CasaXPS `LA(0.9,2,50)`. Satellites are an optional model extension rather than an independently calibrated correction.

Please cite the associated manuscript, *Automated Ce³⁺ Fraction Quantification from Ce 3d XPS Spectra Using Physics-Constrained Bayesian Fitting*, when using this workflow. A final publication DOI can be added when available.

The existing repository [Apache-2.0 license](LICENSE) is retained. Third-party data retain their original terms and attribution; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

The [initial implementation](legacy/v1/README.md) is preserved separately for historical reference. Use the source files at the repository root for the current application.
