# Version 2.0.1 validation — 2026-11-13

This release changes result presentation and file organization. Scientific modules, initialization, proposal construction, convergence diagnostics and sampler defaults are unchanged from version 2.0.0. No long MCMC run was repeated for this update.

- 82 automated tests and 73 subtests passed.
- Both source and frozen Windows applications wrote and read real HTML, PNG, SVG, CSV, JSON and NPZ outputs for LA/Voigt with satellites off/on. Checks covered the primary report, organized directory layout, component areas, posterior quantiles and temporary-file cleanup.
- Report tests covered mean/SD percentage units, retained median/95% interval, section order, selected-scan metadata, missing conditions, HTML escaping, cancelled/failed runs, repeated organization and batch reports.
- GUI tests covered new and legacy result locations, duplicate avoidance, report opening and final-result priority. Hidden geometry checks also covered a completed-run summary at minimum and expanded window sizes.
- A real saved 8000-step result was reformatted without sampling: Ce3+ posterior mean 15.372364%, posterior SD 0.869338 percentage points; the original median and interval were preserved. Its rendered report was visually inspected.
- Full acquisition JSON, complete posterior samples, chain chunks and diagnostics remain available under details/. Unknown acquisition fields are not inferred from survey conditions.

The real-data MCMC and convergence validation of version 2.0.0 remains historical numerical evidence; this update adds reporting/export validation. Quick-fit completion is distinct from passing engineering convergence criteria. Details and executable hashes are recorded in REPORTING_VALIDATION.json. Full local reports are retained in v2.0.1_validation_20261113.
