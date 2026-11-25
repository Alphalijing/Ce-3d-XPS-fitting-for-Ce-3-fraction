# Model definitions and scope

## Component structure and constraints

The ten principal components are ordered `v0, u0, vprime, uprime, v, u, vdouble, udouble, vtriple, utriple`, with nominal centers 880.6, 899.2, 885.2, 903.8, 882.6, 901.2, 888.8, 907.4, 898.4 and 917.0 eV. A common shift is fitted. Spin–orbit partners are separated by 18.6 eV with 3:2 lower-/higher-BE area branching. The Ce(III) pair weights are 2.9:5.4, and the Ce(IV) pair weights are 3.9:2.5:3.8. Source references and physical interpretation are described in the manuscript.

All principal components share Gaussian standard deviation σ. Five kernel-width scales are shared by corresponding spin–orbit partners. Principal hard supports are δE in (-1, 1) eV, σ in (0.05, 1.5) eV, fraction in (0, 1), and kernel widths in (0.1, 2) eV. Full implemented priors and noise terms are tabulated in `data/validation/sampling/default_priors.csv` and defined in `ce3d_xps_model.py`. Width/shift regularization does not turn these supports into universal physical limits for every cerium compound.

Each profile has unit integrated area over the observed fitting window. The Ce³⁺ fraction is the total Ce(III) component area divided by the total Ce component area in that window. With satellites enabled, their Ce(IV) area is included in the denominator.

## Profiles

Voigt profiles are Gaussian–Lorentzian convolutions. Lorentzian γ is the HWHM and Gaussian σ is the standard deviation.

For the LA-type components, the kernel is proportional to `[1 + (x/γ)^2]^(-q)`, with `q = β` for `x < 0` and `q = α` for `x >= 0`, where `x = E - center`. This kernel is convolved with the Gaussian. α = 0.9 and β = 2 are fixed; Gaussian and kernel scales are fitted independently. Only the Ce(IV) v/u components are changed to LA; the other eight remain Voigt. The LA γ is a kernel scale, not the final profile HWHM. This implementation has **no asserted numerical mapping to the CasaXPS third parameter**.

The archived DS comparator similarly modifies only Ce(IV) v/u. Its DS kernel is the expression given in the manuscript Methods, convolved with a Gaussian. Its asymmetry parameter is bounded to [0, 0.05]; finite-window normalization is retained. LA/DS convolution steps are 0.025 eV, with 8192/16384 FFT points, respectively. DS is a restricted comparator, not an optimization over all possible DS extensions, and is not a production GUI option.

## Background and the diagnostic region

The production workflow constructs a noncrossing piecewise Shirley background from the scan ends, adding anchors at crossings, using mean window 1; it then fixes that background. The additional converged-anchor Shirley construction is an exploratory background variation, not an independently measured true background.

The archived regional criterion uses the contiguous low-BE interval within 870–895 eV where the absolute difference between noncrossing and converged-anchor Shirley backgrounds is at most a baseline noise scale. The baseline scale is the median absolute deviation of adjacent intensity differences over 870–875 eV, divided by `0.67448975 * sqrt(2)`. Exact masks are retained with each saved curve. Linear extrapolations are fitted to the measured low-BE baseline. Regional RMSE is an ordinary unweighted residual summary; no regional weight is added to the likelihood.

This diagnostic provides a necessary consistency check within the stated analysis framework, rather than an absolute composition reference. Its interpretation remains conditional on the background variation and spectral coverage considered. Catalytic trends and structural trends provide indirect comparisons, not a known surface Ce³⁺ ground truth.

The survey-derived U4 construction is described as a background with a Tougaard mathematical form. Its fitted lower envelope and intensity transfer do not independently determine a material-specific loss function. Saved curves and fitted parameters are included for inspection.

## Sampling and uncertainty

The likelihood is heteroscedastic Gaussian with `variance = exp(log_a) + exp(log_b) * predicted_peak_signal`. The post-pilot proposal is a fixed mixture of differential-evolution, snooker and a pilot-trained independence Metropolis–Hastings move. Their settings and convergence checks are retained in the sampler configuration. Quick-mode completion is distinct from passing the convergence-validation criteria.

Software headlines use the posterior mean ± 1 posterior SD. Manuscript figures use the posterior median with 2.5% and 97.5% quantiles. These are deliberately different summaries and are identified separately. Credible intervals are conditional on model and background; acquisition differences, beam-induced changes, imperfect templates and background ambiguity can exceed numerical Monte Carlo errors.

## Model extensions and limits

The optional satellite model adds two nonnegative Ce(IV) components constrained to 897–898 and 914.5–915 eV. They share the common Gaussian σ, have fixed Lorentzian HWHM 0.8 eV, and use a joint intensity allocation within the Ce(IV) total. Their centers are on the supplied absolute BE axis. Four additional parameters bring the model from 11 to 15 dimensions. The archived comparison shows that the change can improve some residual regions while worsening others and shifting Ce³⁺; an absolute accuracy gain has not been established.

Joint-background optimization is retained as validation evidence and not exposed as a final production option. Saved paired results show frequent anchor-boundary solutions; one DS joint fit was unsuccessful. These cases are retained and labelled, not silently discarded from the archive.

Strong changes in electronic structure can invalidate fixed multiplet intensity ratios. Unusual sample chemistry, referencing errors, inadequate scan coverage and background/profile mismatch require expert inspection. The included data do not establish a universal calibration for every Ce-containing material.
