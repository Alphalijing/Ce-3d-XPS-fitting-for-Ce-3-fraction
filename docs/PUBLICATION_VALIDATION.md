# Publication-bundle verification — 25 November 2026

The prepared repository passed **82 automated tests and 73 subtests**. Tests were run from the distributed source directory and used the bundled real regression spectrum. A deprecation warning originated in an installed pytest plugin; there were no application test failures.

All ten top-level production Python files are byte-identical to the existing English v2.0.1 source release. The Windows archive matches the recorded SHA-256 and passes ZIP integrity checking. No application rebuilding or long-chain MCMC was performed for this packaging task.

Independent checks verified:

- Six original La-series inputs and backgrounds against the published fit curves.
- All 12 main LA/Voigt Ce³⁺ posterior distributions: 3,840,000 retained draws reproduce the reported medians and 95% bounds.
- The selected 19-spectrum membership, distributed raw hashes and original-source lineage, saved backgrounds and contiguous-region residuals.
- All 48 synthetic paired inputs and known generating backgrounds.
- Regeneration of all 48 synthetic cases: **144 raw, truth-curve and component CSV files are byte-identical** to their archives.
- The independent-parameter LA prediction and analytic Jacobian against numerical checks.
- Recalculation of the reported model-comparison, synthetic-error and joint-background statistics.
- Successful execution of the Figure 3, 4, 5 and six-row posterior plotting scripts using only bundled data.

The recreated figures were inspected as images. File hashes and relative documentation links were checked for the final package. Machine-readable results and dependency versions are recorded in `publication_validation.json`.

The original reporting/export checks remain described in the root `VALIDATION.md`. Model plausibility and absolute chemical accuracy are separate from these numerical and software checks; the associated conditional scope is described in `MODEL_AND_SCOPE.md`.

Before public upload, 26 XLS copies received an equal-length redaction of an identifying operator account inside embedded Windows paths and workbook author/write-access metadata. Every worksheet cell was checked against the source; spectral numbers and acquisition conditions are unchanged. Both byte hashes are retained in `workbook_redaction.json`.
