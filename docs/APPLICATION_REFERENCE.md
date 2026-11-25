# Ce 3d XPS Fitter 2.0.1 — English Edition

Physically constrained Bayesian fitting for Ce 3d XPS. The default model uses Gaussian-convolved LA-type peaks (fixed alpha = 0.9, beta = 2) and a fixed Shirley background. Voigt peaks and Ce(IV) satellite components are optional.

Version 2.0.1 puts the Ce3+ result first and groups detailed records into subfolders. The scientific model, numerical algorithms and sampling defaults remain unchanged from version 2.0.0.

## Start the application

**Windows package:** extract the complete archive and double-click `Ce3d_XPS_Fitter_v2.0.1_EN.exe`. Keep the adjacent `_internal` folder; the executable cannot be distributed alone.

**Python source:** Python 3.11 is recommended. Run `python -m pip install -r requirements.txt`, then double-click `Launch_GUI.bat` or run `python ce3d_xps_fitter_gui.py`. See [GUI_GUIDE.md](../GUI_GUIDE.md) for the interface guide.

## Defaults and sampling modes

| Setting | Default |
|---|---|
| Peak shape | LA; fixed alpha = 0.9, beta = 2; fitted Gaussian width |
| Background | Shirley background constructed first, then fixed |
| Satellites | Off |
| Mode | Quick fit |
| Walkers | 80 |
| Total steps | 8000, including a 2000-step pilot |
| Discard | First 4000 steps, including the pilot |
| Progress messages | At start, every 2000 steps, and on completion or stop |

**Quick fit** runs MCMC. Least-squares fitting and multistart posterior optimization provide initial positions, followed by the validated sampling strategy. It is not a least-squares-only mode.

**Convergence validation** checks every 2000 steps after burn-in, with no fixed upper step limit by default. Convergence cannot be declared before 10000 total steps. Sampling continues until all engineering criteria pass or the user stops it:

- Retained length is at least 100 times the slowest integrated autocorrelation time, including that of the log posterior.
- Across the latest three checkpoints, both adjacent relative changes in every autocorrelation time are below 5%.
- Across the latest three non-overlapping 2000-step windows, the ranges of the Ce3+ median and each 95% interval endpoint are each below 0.1 percentage points.

**Stop** saves completed chain blocks and available statistics after the current sampling segment finishes. Cancellation, quick completion and validated convergence have distinct statuses. These checks assess numerical sampling under the selected conditional model; they do not establish that the background or peak shape is uniquely correct.

## Input

- **`.xls/.xlsx`:** automatically locates `Ce3d Scan`, accepting common spacing and case variants. Column A is binding energy and column C is raw intensity; a survey sheet is optional. Explicitly choosing the provided background uses column E.
- **`.csv/.txt/.dat/.xy`:** two columns for binding energy and raw intensity; common comma, tab and whitespace delimiters are supported. Provided-background mode uses column three. Header rows are allowed; missing values and duplicate energies in the numeric data are rejected.
- **`.vms/.vamas`:** supports REGULAR XPS energy scans. A unique Ce 3d block can be detected automatically; multiple Ce 3d blocks or signal channels require an explicit selection. Survey and transmission channels are not treated as Ce 3d intensity.
- **VAMAS kinetic-energy axes:** require an explicit offset for `BE = offset - KE`. The program does not infer vendor work-function or energy-reference conventions. This offset is unused for binding-energy input.
- By default, the input's low-energy end is retained and the upper end is cropped to 925 eV. The validated adaptive, non-crossing Shirley procedure starts at both ends with mean window = 1, then fixes the resulting background. Input must fully cover 879–918 eV; a broad survey cannot replace a narrow scan.

Instrument intensity values are used as supplied. Transmission correction, charge referencing and cross-scan intensity calibration are not applied automatically. Available experimental conditions and explicit input transformations are recorded in JSON. Typical inputs use counts or counts/s; posterior area and noise parameters inherit the corresponding intensity units.

## Experimental metadata

Each run writes `details/experiment_metadata.json`. The main report shows a concise summary of the selected Ce 3d acquisition conditions. Excel metadata includes available acquisition time (for example, cell I4), pass energy, energy step and original cell locations. Ce 3d and survey conditions are stored separately. VAMAS metadata preserves the header, block inventory, selected block/channel, acquisition conditions, energy-axis transformation and source SHA-256. Missing or unspecified conditions remain `null`; file format alone does not establish transmission correction or charge calibration.

## Peak shapes, constraints and area definition

The ten main components retain constrained relative positions, relative area weights within Ce(III)/Ce(IV), a 3:2 doublet branching ratio and an 18.6 eV separation. Components share Gaussian sigma; paired components share a width parameter gamma. For Voigt, gamma is Lorentzian HWHM. For LA, gamma is the base function's scale parameter, not the actual HWHM of an asymmetric or convolved peak; output tables identify the distinction.

In LA mode, only Ce(IV) v/u use LA peaks; the remaining main components use Voigt peaks. Voigt mode uses symmetric Voigt peaks throughout. Every component is normalized by its integral within the measured energy window. The fixed-alpha/beta Gaussian-convolved LA implementation does not claim to reproduce a particular CasaXPS third parameter.

The heteroscedastic Gaussian likelihood uses `variance = exp(log_a) + exp(log_b) * peak_model`. The main model retains the validated priors and parameter bounds. The former regional confidence-weight parameter has been removed. Background stays fixed during MCMC; unstable joint background optimization is not part of the default workflow.

**Include satellites** adds two candidate Ce(IV) Voigt components. Their absolute centers are constrained to 897–898 and 914.5–915 eV, without the main peaks' global dE shift. They share Gaussian sigma; their Lorentzian HWHM defaults to a fixed 0.8 eV. Their nonnegative intensities are allocated jointly within the total Ce(IV) area, without an intensity ratio unsupported by a reference sample. These settings and constraints are recorded in JSON.

With satellites enabled, the Ce3+ fraction is the measured-window Ce(III) area divided by the measured-window area of all Ce components, including Ce(IV) satellites. This option supports model-sensitivity comparisons; it does not imply that every sample requires resolvable satellites or that adding them is always more accurate.

Satellites increase the model from 11 to 15 dimensions. The main 8000-step precision validation concerns the default model without satellites. Use Convergence validation to examine the more complex satellite posterior further.

## Results first

Click **Open report** after a run, or open `Results.html` in its results folder. The report is ordered for review:

1. **Ce3+ percentage ± uncertainty**, displayed as the posterior mean ± one posterior standard deviation. The uncertainty is in percentage points. The posterior median and exact 95% credible interval are shown separately; ±1 SD is not the 95% interval.
2. **Component fit plot**, followed by a concise summary of the selected Ce 3d scan's recorded acquisition conditions. Missing conditions are not inferred from the survey scan.
3. **Per-component fitting results**, including centers, Gaussian sigma, correctly identified Voigt/LA width parameters and measured-window areas.
4. **Ce3+ posterior distribution**, with the remaining parameter distributions and full posterior statistics available below.

The report embeds its figures and can be copied as a single HTML file. The accompanying folders retain exportable figures, tables and complete records:

| File or directory | Contents |
|---|---|
| `Results.html` | Main report; a batch report lists all sample results and links to their individual reports |
| `figures/` | Component/background plots and posterior distribution figures, with PNG/SVG exports |
| `tables/` | `peak_table.csv`, `fit_curves.csv` and `posterior_parameters.csv` |
| `details/` | Full result JSON, experimental metadata, configuration, retained posterior samples, prepared spectrum, initialization, checkpoints, progress and closed chain blocks |

The interface's final result line also shows the primary Ce3+ estimate before output paths. Full experimental metadata remains in `details/experiment_metadata.json`; numerical summaries remain in `details/fit_result.json`. Detailed records are retained for reproducibility without crowding the first view.

Finishing 8000 quick-mode steps does not imply that all engineering convergence criteria passed. Posterior uncertainty is conditional on the selected fixed background, peak shapes, constraints, priors and noise model; it does not include all experimental and model systematic errors.

## Command-line examples

```powershell
python ce3d_xps_fitter_core.py --input examples/La10at_validation.xls --outdir results
python ce3d_xps_fitter_core.py --input examples/La10at_validation.xls --outdir results --model Voigt
python ce3d_xps_fitter_core.py --input examples/La10at_validation.xls --outdir results --satellites
python ce3d_xps_fitter_core.py --input examples/La10at_validation.xls --outdir results --mode convergence
python ce3d_xps_fitter_core.py --input examples/GDC10_published.vms --outdir results --vamas-block 12 --vamas-ordinate 0 --kinetic-offset-ev 1486.6
```

The final example's block, channel and offset apply only to the bundled public GDC example. See [examples/README_input_sources.md](../examples/README_input_sources.md) for provenance and source links. `La10at_validation.xls` is the project's La-labelled 0.5 sample, corresponding to 10 at% doping. `validation_sample_0p5.csv` contains its real data and validation background for model regression; neither is a synthetic spectrum or an independent physical standard.

## Tests and Windows build

```powershell
python -m pip install -r requirements-dev.txt
python -m pytest tests -q
python build_windows_exe.py
```

The build creates `dist/Ce3d_XPS_Fitter_v2.0.1_EN/Ce3d_XPS_Fitter_v2.0.1_EN.exe`. Distribute the entire application folder.

See [VALIDATION.md](../VALIDATION.md) for the inherited scientific validation and separate English-edition checks. Model tests compare the validated LA/Voigt posteriors element by element and cover satellite areas, mixed-MH correction, chain continuity, input parsing, run status and cancellation. Version 2.0.1 does not modify original input files.
