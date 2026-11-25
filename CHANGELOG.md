# 2.0.1 — 2026-11-13

## Results-first report

- Added a self-contained `Results.html` with the Ce3+ posterior mean ±1 SD first, followed by the fit figure, concise selected-scan conditions, component parameters and posterior distributions. Median and 95% credible interval are retained separately.
- Grouped figure exports, tables and complete technical records under `figures/`, `tables/` and `details/`.
- Changed Open results to Open report and placed the final Ce3+ summary before output paths in the interface log.
- Added a batch report linking individual sample reports. Preserved legacy status-reading support.
- Kept the scientific model, sampler, defaults, component-area definitions and saved-chain values unchanged.

# 2.0.0 — 2026-11-12

## SVG export packaging fix — 2026-11-12

- Included the SVG backend explicitly and verified its availability before sampling.
- Added full PNG/SVG/CSV/JSON export self-checks to the frozen executable, covering LA/Voigt with satellites off/on without rerunning MCMC.
- Fitting models, priors, sampler settings and saved-chain values are unchanged.

## English edition

- Translated the interface, launcher messages and user documentation into English.
- Added `_EN` to the Windows application folder and executable names.
- Retained version 2.0.0, scientific algorithms, model constraints, defaults and input/output behavior. This is a language edition, not a new scientific model release.
- Inherited the original scientific validation; English-edition checks are recorded separately in `VALIDATION.md`.

## Original 2.0.0 changes

- Made LA peaks with a fixed Shirley background the default; kept Voigt as an option.
- Set defaults to 80 walkers, 8000 total steps and 4000 discarded steps. Mixed-MH proposals are frozen after the 2000-step pilot and combined with DE/Snooker sampling.
- Added Convergence validation with no fixed step limit by default, plus Stop and save.
- Added optional Ce(IV) satellites, off by default, with recorded center, width and intensity constraints and a satellite-inclusive Ce3+ area definition.
- Removed regional confidence weights and their interface controls.
- Added Excel, common text and VAMAS input detection, with explicit VAMAS block/channel selection.
- Preserved available Excel/VAMAS acquisition metadata, provenance and unknown values in JSON.
- Arranged model, experimental conditions and run mode in three option rows, followed by input, output and the run log; added start and 2000-step progress feedback.
- Added chain-block saving, parameter/component exports, distinct run statuses, cancellation and batch status.
- Added real examples, regression/sampling tests and a Windows one-directory application.
