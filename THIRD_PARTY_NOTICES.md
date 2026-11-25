# Third-party data and dependencies

## Public VAMAS example

`examples/GDC10_published.vms` is an unchanged copy of the publicly supplied supplementary file `ce0.9gd0.1o2-d _gdc10_ xps.vms` from:

A. Nenning, S. Breitwieser, C. Melcher and J. Fleig, *UHV-based analytics with electrochemical oxygen activity control*, Journal of Materials Chemistry A **13** (2025) 29147–29160. [doi:10.1039/D5TA02648B](https://doi.org/10.1039/D5TA02648B).

The [article and its supplementary-material links](https://pmc.ncbi.nlm.nih.gov/articles/PMC12326307/) identify a Creative Commons Attribution 3.0 Unported licence. The article is © The Royal Society of Chemistry. Retain this attribution when redistributing the example; the repository's Apache licence does not replace its original terms. [Original supplementary archive](https://pmc-oa-opendata.s3.amazonaws.com/PMC12326307.1/TA-013-D5TA02648B-s002.zip).

The file contains Gd-doped ceria under electrochemical conditions. It is an input-compatibility example, not a certified pure Ce(IV) standard. Its SHA-256 is `57377a8f817613bbea1705e13b50abb20bb5597b873ac3e74f5dd940a32d2111`.

## Dependencies

NumPy, SciPy, pandas, Matplotlib, emcee, vamas, openpyxl, xlrd and tkinterdnd2 are separate projects under their respective licences. Python packages are installed using the requirements files. The Windows distribution retains the accompanying dependency files from the validated build.

No third-party journal PDFs or proprietary CasaXPS software are included.
