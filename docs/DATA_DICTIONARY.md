# Data dictionary

All binding energies and widths are in eV. Intensities retain the source instrument scale (usually counts s⁻¹ for the experimental series). A percentage-point change is a difference between percentage values, rather than a relative percent change.

| File group | Principal fields |
|---|---|
| `la_series/sample_index.csv` | Original label, publication sample ID, nominal La at%, relative raw/prepared paths and file SHA-256 |
| `la_series/curves/*.csv` | `energy_eV`, `raw`, `high_mask`, and model-specific total/background/residual curves; residual = raw − total |
| `la_series/Figure4_all_model_results.csv` | Model and estimator; Ce³⁺ median/interval where available; shared σ and width quantiles; total and regional RMSE |
| `la_series/Figure4_composition_panels_data.csv` | Six-row table for model residuals, Ce³⁺, Gaussian σ and DMC productivity used in Figure 4 |
| `la_series/posterior/*.npz` | `Ce3_fraction`: all post-discard Ce³⁺ draws (fraction units); `parameter_draws_temporal_subset`: visualization subset; parameter names and retained step range |
| `additional_19/summary.csv` | Source membership; raw-file hash; archived spectral-content hash; contiguous interval endpoints/count; LA/Voigt regional RMSE and relative reduction; boundary flags |
| `additional_19/curves/*.csv` | Raw data, noncrossing/converged-anchor Shirley constructions, exact contiguous high-confidence mask and archived fit curves |
| `experimental/*/prepared/*.csv` | Three columns: ascending energy, raw intensity, saved fixed background; select provided-background mode to retain column three |
| `synthetic/truth/*/*.csv` | Descending-energy simulated observations with their true background; associated `_truth_curves` and `_components` files use ascending energy |
| `synthetic/truth/*/*_truth.json` | True composition, widths, component areas, individual Poisson seed and generation details |
| `synthetic/prepared/*.csv` | Same known background for both fit methods; cropped 545-point fitting window |
| `synthetic/LA_MCMC_summary.csv` | Posterior median, mean, SD, q025/q975, signed errors in percentage points, RMSE and sampling status |
| `synthetic/free_LA_summary.csv` | Independent-parameter least-squares estimate, RMSE, chosen start and constrained-LS diagnostic seed result |
| `synthetic/fits/` | Individual fitted curves/peak tables, parameter summaries and run configuration; large full chains are omitted |
| `validation/background/` | Six saved background variations; exact regional masks; local linear extrapolation; U4 parameter table and transferred curves |
| `validation/satellites/` | Paired curves, Ce³⁺/residual comparisons, constraints, posterior parameters and original sampling status |
| `validation/joint_background/` | Paired deterministic fixed/joint results, optimizer settings, bound flags and saved backgrounds/curves |
| `validation/sampling/` | Prior table, sampler settings, retained-prefix quantiles, autocorrelation/precision histories and seed comparison |
| `validation/prior_sensitivity/` | Soft-prior variations and a separately labelled hard-support expansion; shifts are conditional comparisons |

`sigma` denotes Gaussian standard deviation, not Gaussian FWHM. For a Voigt component, `gamma` is Lorentzian HWHM. For LA/DS components it is a kernel scale, not the final asymmetric peak's HWHM. `Ce3_low/high` and `q025/q975` refer to 95% equal-tail credible bounds where posterior sampling is available. DS point estimates have no artificial interval added.

Some source tables retain historical run identifiers to trace paired analyses. Internal algorithm labels such as `v7` are preserved only in original archived evidence where necessary; the manuscript-facing descriptions use noncrossing/fixed Shirley and converged-anchor Shirley. Old evidence-file numbering is replaced by descriptive directory names.

Numerical configuration, source-byte hashes and third-party attribution are provided separately. Missing instrument fields are not inferred.
