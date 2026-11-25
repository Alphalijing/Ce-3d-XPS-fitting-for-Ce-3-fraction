# Publication release bundle — 25 November 2026

The application version remains **2.0.1**. Production scientific modules and GUI source are copied without changes from the validated English release. The Windows executable and dependencies are unchanged; one bundled example workbook has identifying operator metadata redacted, matching the source/data package. This bundle adds publication data, portable analysis scripts, source/data hashes and updated documentation.

Relative to the initial manuscript's software, the workflow uses LA-type peaks and a fixed noncrossing Shirley background by default; Voigt and satellite options are available; the former regional confidence parameter is removed; quick MCMC defaults to 80 walkers and 8000 steps; an optional convergence-validation mode is supplied. Input handling includes native VAMAS scans and experimental-condition JSON output. Version 2.0.1 presents Ce³⁺ and uncertainty first and groups detailed outputs.

The regression test's one personal-machine fallback path has been removed from this upload copy. Its bundled real validation spectrum and frozen reference calculations are retained. No scientific parameter or fitting result was changed for packaging.

The numerical results retain their original estimators and backgrounds. The archived extended screen, synthetic known-background benchmark and satellite/joint paired comparisons are labelled separately from the main experimental MCMC comparison.
