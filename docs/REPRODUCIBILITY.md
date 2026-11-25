# Reproducing the publication results

All paths below are relative to the repository root. Python 3.11 and the supplied requirements are recommended. Hashes in `MANIFEST.json` identify the distributed files. Spectral values and acquisition conditions are preserved. Public Excel copies have an identifying operator account redacted from embedded Windows source paths and workbook author/write-access metadata; this equal-length metadata substitution is documented in `workbook_redaction.json`. The original and distributed workbook hashes are both retained. The already-published VAMAS example remains byte-identical to its source. Exported metadata and evidence tables use portable paths or historical run identifiers rather than personal-machine paths.

## What is included

| Manuscript location | Distributed evidence | Reproduction |
|---|---|---|
| Section 3.1 / Figure 3 | `data/validation/background/La_2at.csv`; six background constructions; high-confidence-region table | `python analysis/plot_figure3.py` |
| Section 3.2 / Figure 4 | `data/experimental/la_series/`; six raw workbooks, exact masks/curves, composition/catalysis data and posterior draws | `python analysis/plot_figure4.py` |
| Section 3.2 extension | `data/experimental/additional_19/summary.csv`, raw/prepared inputs and saved fits | `python analysis/summarize_results.py` |
| Section 3.3 / Figure 5 | `data/synthetic/`; 48 generated and prepared inputs, truth, both fitted workflows and 95% intervals | `python analysis/plot_figure5.py` |
| SI posterior figure | `data/experimental/la_series/posterior/LA_La_*at.npz` | `python analysis/plot_posterior_s2.py` |
| SI background, satellites and joint backgrounds | `data/validation/background`, `satellites`, `joint_background` | Inspect supplied curves/configuration and summary tables |
| SI sampling and prior sensitivity | `data/validation/sampling`, `prior_sensitivity` | Inspect saved checkpoint/precision tables and prior variants |
| SI input examples | `examples/` and `data/validation/input_examples/` | See commands below |

The original La-series labels `0, 0.05, 0.1, 0.2, 0.5, 1` correspond to **0, 1, 2, 4, 10, 20 at% La**. Thus the original label `0.5` denotes 10 at% La. The mapping is explicit in `sample_index.csv`.

The additional dataset contains `batch2_1`, `batch2_3`–`batch2_7`, and `batch2_9`–`batch2_21`. It excludes the separate F-series and two repeated La-series spectra. These are 19 additional spectra, without a claim that all represent independently prepared specimens. Instrument conditions remain in the raw files and extracted metadata; unspecified sample identities are not inferred.

## Existing estimates versus new fits

LA and Voigt main-text results use 80 walkers, 8000 total steps and a 4000-step discard. DS remains an archived least-squares point estimate, with no MCMC interval assigned to it. Regional RMSE is computed from the saved fit curve at the same measured points and fixed background for each specimen. It is not a likelihood-weighted metric.

The `.npz` posterior files contain **all 320,000 retained Ce³⁺ draws** used for each of the 12 LA/Voigt main estimates. `Ce3_fraction` is a fraction, not percent. Each file additionally contains a temporal subset of the 11-dimensional parameter draws (all 80 walkers every 100 retained steps) for visualization. The latter subset is not used to recompute the tabulated full-chain parameter quantiles. Full 11-dimensional chains and full log-probability arrays are not duplicated in this compact bundle; the reported parameter quantiles and original source-chunk hashes are retained.

The extended 19-spectrum screen used archived deterministic fits and backgrounds produced before the production preprocessor was frozen. Reconstructing a background with today's default can differ from those saved backgrounds. **Use `additional_19/prepared/*.csv` with `--background-mode provided` to retain the exact background of that screen.** The size of each reconstruction difference is recorded in the summary and metadata. Neither the stored fit curves nor the published residual metric has been replaced by a newly computed result.

The spectral-content hash in `additional_19/summary.csv` is distinct from the raw-workbook byte hash; both are retained and labelled.

Example main-series fit using the released default:

```console
python ce3d_xps_fitter_core.py --input data/experimental/la_series/raw/La_10at.xls --outdir results/new_La10_fit
```

To keep an archived background exactly:

```console
python ce3d_xps_fitter_core.py --input data/experimental/additional_19/prepared/batch2_1.csv --background-mode provided --outdir results/new_extension_fit
```

New MCMC runs are stochastic and need not reproduce the last reported digits. The plotting and verification commands only use the saved results.

## Synthetic data

The 48 cases cover Ce³⁺ fractions of 0, 15, 30 and 60%, Gaussian σ of 0.5, 1.0 and 1.5 eV, and nominal Poisson photo-peak noise settings of 3 and 10%, in two background conditions. The generating photoelectron envelope uses the same fixed-α/β LA-type model as the constrained fit. This is a within-model recovery benchmark.

On the increasing binding-energy axis, define `Q(E) = cumulative_trapezoid(P(E), E, initial=0)`. The added cumulative tail is `S(E) = max(P) Q(E)/max(Q)`. The second condition uses `B(E) = 0.86 max(P) + 0.84 S(E)`; the first uses `B = 0`. Poisson draws use expectation `P + B`. No decay kernel and no clipping of observations to the background are applied. Full generation parameters, individual seeds, truth curves and component curves are supplied.

**Both fitted workflows use the same known generating background**, fixed in `data/synthetic/prepared/`. They do not reconstruct it with Shirley. The constrained workflow uses MCMC; the comparator uses 40 independent component parameters with ordinary least squares and four starts. Comparison therefore concerns the complete constrained Bayesian workflow against that independent-parameter fit, rather than the sampler alone.

Regenerate the raw/truth curves without fitting:

```console
python analysis/generate_synthetic.py
```

The generated files are written into `results/regenerated_synthetic/`; distributed truth data are preserved. Random seeds and CSV numeric formatting are retained. The generation grids extend to 930.58 eV; fitting crops to 924.98 eV, retaining 545 points. The component area convention refers to the fitting window.

Optional refitting of one synthetic case (a full MCMC fit takes time):

```console
python analysis/refit_synthetic.py --case a_LA_Ce3_015pct_sigma_1eV_poisson_03pct --method MCMC
python analysis/refit_synthetic.py --case a_LA_Ce3_015pct_sigma_1eV_poisson_03pct --method free_LA
```

The wrapper uses the recorded case seed and fixed supplied background. All new results go to `results/synthetic_refits/`. `python analysis/free_LA_fit.py --check` checks the independent-parameter LA model and its analytic Jacobian without fitting cases.

## VAMAS input compatibility

The attributed public GDC example includes 61 blocks and multiple channels. For the example used here, the selected Ce 3d block is 12 (zero-based), the ordinate is 0, and the explicitly chosen transformation is `BE = 1486.6 - KE`:

```console
python ce3d_xps_fitter_core.py --input examples/GDC10_published.vms --vamas-block 12 --vamas-ordinate 0 --kinetic-offset-ev 1486.6 --outdir results/GDC_example
```

This command now uses the release defaults (LA, 8000 steps). The archived operational smoke test used **Voigt, 4000 total steps and 2000 discarded steps**. Its result is an import/export check and is not a certified Ce³⁺ composition or a convergence claim. The selection and offset are example-specific and must not be applied automatically to other VAMAS files.

## Further checks

`python analysis/verify_release.py` checks file hashes, all main Ce³⁺ posterior quantiles, regional RMSE, selected dataset membership, paired synthetic inputs and key statistics without sampling. Its JSON output and `python analysis/summarize_results.py` output go to `results/`.

Satellite and joint-background evidence retains its original paired baselines. These are not substituted for the separate main-text MCMC series. Satellite parameters add four dimensions and are not fixed intensity ratios calibrated against a pure-Ce(IV) standard. The joint-background comparison retains its optimization status and boundary-hit flags, including the unsuccessful DS case.

## Catalytic-data provenance

The DMC-productivity values in Figure 4 were transcribed from the numerical labels printed above the bars in Figure 6 of our previous study: L. Hu, S. Zhang, P. Zhao, J. Yang, G. Wu and W. Xu, *La-doped CeO₂ nanoflowers catalysts for direct synthesis of dimethyl carbonate from CO₂ and methanol*, Journal of Rare Earths **43** (2025) 2177–2185, [doi:10.1016/j.jre.2025.03.026](https://doi.org/10.1016/j.jre.2025.03.026). The six values are 7.25, 8.70, 8.81, 8.16, 9.42 and 2.05 mmol g⁻¹, in increasing La composition. These are printed numerical values rather than coordinates estimated from a plotted curve. The transcribed table and source provenance are retained in `data/experimental/la_series/catalysis/`. The original journal PDF is not redistributed.

## Public-workbook privacy preparation

Twenty-six distributed XLS files (25 experimental workbooks and the repeated La10at example) contain a privacy-only substitution in embedded source-file paths and workbook author/write-access metadata. Every worksheet cell was checked against the original workbook: all numbers, spectra and acquisition conditions remain unchanged; identifying path strings, where visible as cells, use a generic account. File lengths and workbook structure are preserved. `docs/workbook_redaction.json` lists original-source and public-copy SHA-256 values. `raw_file_sha256` in the dataset indexes identifies the distributed copy; `original_workbook_sha256` records its pre-redaction source. Historical fit-run hashes refer to those original sources. The original local workbooks were not modified.
