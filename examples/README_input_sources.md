# Raw input examples and provenance

## La10at validation data

`La10at_validation.xls` is the project's real La-labelled 0.5 experimental sample, corresponding to 10 at% doping. It is the raw Excel input for the reference LA runs in `REFERENCE_RESULTS.json`.

SHA-256: `10b5a67b602faa006532e763da98cc9041d98e43679b2da4c547971f58acc6c8`.

`validation_sample_0p5.csv` contains data from the same sample and its validation background for model-regression tests. It is not a synthetic spectrum, a separate experimental sample or an independent certified valence standard.

SHA-256: `61845818b9d81828be40ead28402f894041c816be8bcd25cccaf6eb4cf3ebb52`.

## Commercial CeO2 Excel example

`shang_yong_CeO2_raw.xls` is an unchanged copy of the user's experimental file `shang yong-CeO2.xls`.

SHA-256: `654cb1786d0e8504172b630d4cacf0eba99e9ecb400c0ed3f9b1e20151d58644`.

It contains `Ce3d Scan`, other scans, a survey and vendor-export tables. Ce 3d acquisition metadata is read from the Ce 3d sheet only. Cells H4/I4 explicitly give `Total acquisition time` / `8 mins 0.4 secs` (480.4 s). Ce 3d pass energy is 50 eV; survey pass energy is 150 eV. The export does not explicitly give a numeric photon energy, transmission-correction status or charge reference, so these fields remain `null`. The vendor background in column E is used only when provided-background mode is explicitly selected.

## Public GDC10 VAMAS example

`GDC10_published.vms` is the unchanged public experimental file `ce0.9gd0.1o2-d _gdc10_ xps.vms`, retained from the project's earlier Ce VAMAS import test.

SHA-256: `57377a8f817613bbea1705e13b50abb20bb5597b873ac3e74f5dd940a32d2111`.

Source: Nenning, Breitwieser, Melcher and Fleig (2025), *UHV-based analytics with electrochemical oxygen activity control*, [DOI: 10.1039/D5TA02648B](https://doi.org/10.1039/D5TA02648B). The experimental file is available in the [public supplementary ZIP (PMC Open Access mirror)](https://pmc-oa-opendata.s3.amazonaws.com/PMC12326307.1/TA-013-D5TA02648B-s002.zip); consult the source publication and supplementary material for reuse terms.

This is Gd-doped ceria, not an undoped certified CeO2 valence standard. Its 61 blocks include 11 Ce 3d scans and multiple ordinate channels. The bundled import example uses block 12 (`+800 mV`), ordinate 0 (`Intensity`) and the explicitly supplied conversion `BE = 1486.6 - KE`, giving 874–924 eV. These values are not selected automatically. Confirm the instrument's energy convention before using an offset in scientific analysis.

The transmission ordinate is neither fitted as intensity nor silently applied as a correction. Existing CasaXPS comments are preserved as source metadata and are not executed.

## Scope of the examples

Tests also construct clearly named temporary synthetic text/Excel inputs and a synthetic axis-label edit to check validation. These are parser tests, not new experimental samples or MCMC results.

Import examples exercise file parsing and background preparation. The reference La runs additionally test fitting as described in `../VALIDATION.md`; bundled examples alone do not establish model validity, certified composition or convergence for other data. The short public-VAMAS run is an integration check with a custom budget.
