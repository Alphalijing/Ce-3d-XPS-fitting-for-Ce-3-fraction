# Automated Ce 3d XPS Fitter

This tool performs automated, physically constrained fitting of Ce 3d XPS spectra and reports the Ce3+ fraction with MCMC-derived uncertainty. It is designed for Ce 3d analysis only. O 1s fitting is not included in this release.

## Main Features

- Windows GUI for routine use.
- Command-line interface for scripted workflows.
- Raw XPS Excel input mode.
- Processed CSV/table input mode.
- Physics-constrained ten-component Ce 3d model.
- Shirley-background-derived confidence weighting for the Gaussian likelihood.
- MCMC posterior sampling using `emcee`.
- Outputs fit curves, peak tables, posterior samples, parameter summaries, and a JSON summary.

## Input Mode 1: Raw XPS Workbook

Use this mode for the original XPS workbook exported by the instrument/software.

Required sheet names:

- `XPS Survey`
- `Ce3d Scan`

The program uses the sheet names to avoid fitting the wrong spectral region.

For `Ce3d Scan`, the program reads:

- Column A: Binding energy
- Column C: Raw intensity
- Column E: Background

The first numeric background value in raw XPS exports is often `0`. The program treats this value as a software placeholder and drops that first data row before background subtraction.

The fitted intensity is:

```text
background-subtracted intensity = raw intensity - background
```

## Input Mode 2: Processed CSV / Table

Use this mode when the Ce 3d spectrum has already been exported to a simple three-column file.

Required column order:

1. Binding energy in eV
2. Raw intensity before background subtraction
3. Background intensity

The program subtracts the third column from the second column before fitting:

```text
background-subtracted intensity = raw intensity - background
```

The processed file should contain one Ce 3d spectrum only. Header rows are tolerated only when the numeric data can still be recognized, but a clean header-free CSV file is recommended. The Ce 3d binding-energy region is typically around 870-930 eV.

## Likelihood Weighting

The Ce 3d likelihood uses a Shirley-background-derived confidence factor:

```text
C(E) = 1 / [1 + lambda * T(E)]
```

where `T(E)` is the three-segment scattering-tail accumulation index derived from the Ce 3d background.

Available GUI options:

- `lambda = 0`: unweighted Gaussian likelihood.
- `lambda = 2`: recommended default.
- Custom non-negative lambda value.

## Default MCMC Settings

- walkers: 80
- steps: 5000
- burn: 1500
- thin: 10
- seed: 42

## Output Files

Each run creates a timestamped folder under the selected output directory. Typical outputs include:

- `fit_result.json`
- `fit_ce3d.png`
- `fit_curve_ce3d.csv`
- `peak_table.csv`
- `peak_table_ce3d.csv`
- `posterior_samples_ce3d.csv`
- `posterior_param_summary_ce3d.csv`
- `prepared_input/processed_ce3d_bg_sub.csv`
- `prepared_input/ce3d_background.csv`
- `prepared_input/ce3d_raw_with_background.csv`

## Running from Python

Install the dependencies:

```bash
pip install -r requirements.txt
```

Run the GUI:

```bash
python ce3d_xps_fitter_gui.py
```

Run from the command line:

```bash
python ce3d_xps_fitter_core.py --input path/to/raw_xps.xls --raw-xps --lambda-ce 2 --outdir output_folder
```

For processed CSV input:

```bash
python ce3d_xps_fitter_core.py --input path/to/ce3d.csv --lambda-ce 2 --outdir output_folder
```

## Model Notes

The Ce 3d envelope is modeled using ten physically constrained Voigt-type components corresponding to Ce3+ and Ce4+ final-state features. The model applies fixed spin-orbit splitting, degeneracy-derived 3d5/2:3d3/2 intensity ratios, constrained Ce3+/Ce4+ multiplet patterns, grouped Lorentzian broadening, and a shared Ce 3d Gaussian broadening parameter.
