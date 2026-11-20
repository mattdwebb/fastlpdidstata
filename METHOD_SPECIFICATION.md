# DGJT LP-DiD Reference Specification

Reference: Dube, Girardi, Jorda and Taylor (2025),
https://doi.org/10.1002/jae.70000, sections 3 and 4.1.
Author examples frozen at 5d7809c8d818ac2ba4d1d53bf38696763a77544f.
Official Stata lpdid 1.1.0 (20sep2026) is the executed external oracle.

For each base time t and horizon h, regress Y(i,t+h)-Y(i,t-1) on
D(i,t)-D(i,t-1) and time effects. Keep new entrants or controls with
D(i,t+max(h,0))=0. Require actual calendar observations for the needed lags,
leads, covariates, treatment difference, and outcomes. Gaps are never compressed.
Treatment must be absorbing, binary, and nonmissing. Cohort zero/missing is never.
Treatment-only adoption in an unobserved gap is rejected; supply the known cohort.

Variance weighting uses OLS. Reweighted (default) uses WLS weights n_t/n_control,t
on the final complete estimation sample, excluding cells with no treated or
no controls. These weights make each contributing treated entry equally weighted.
Variance weighting removes singleton time cells, as reghdfe does.

Pooled scalar uses mean(Y(i,t+a),...,Y(i,t+b))-Y(i,t-1) and controls clean
through b. All outcomes in the pooled window must exist. It is not an average
of separately fitted horizon coefficients. event_window=(lo,hi) is mandatory:
pooled post window max(0,lo)..hi, horizon estimates lo..hi, -1 reference.

With covariates, reweighted uses the control-only regression
Z = time effects + X gamma, then averages Z-X gamma-time effect over entrants.
Retain all-control time cells when estimating slopes, but exclude cells with no
controls. Singular designs fail explicitly. Variance weighting includes additive
covariates in OLS and needs stronger assumptions about effect heterogeneity.
Covariates default to t-1. Explicit lag=0 uses caller-supplied predetermined
entry-row covariates; no automatic freezing and no post-treatment covariates.

## Linearity

Conditional on fixed treatment, sample, covariate matrix, and weights, WLS is
e1'(X'WX)^(-1)X'W L Y, where L is the outcome-difference/pooled transformation.
RA is mean_T(LY) - mean_T(X)(X_C'X_C)^(-1)X_C'L_C Y.
Both are exactly linear in Y. The numerical suite checks superposition.
If covariates are reconstructed from Y (for example lagged outcomes), the
fixed-design claim no longer applies. Assignment changes require rebuilding
clean samples, time cells, weights and all fits; no frozen-design shortcut is used.
