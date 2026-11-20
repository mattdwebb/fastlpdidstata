# Validation Status: Partial

Executed official SSC lpdid 1.1.0 in StataNow 19.5 SE on deterministic synthetic
panels. Native Stata point estimates and VCEs pass at absolute tolerance 1e-9
for variance/reweighted, with/without covariates, pooled and horizons -2,0,1,2.
Base fixture: 20 point + 20 SE checks. Missing outcome and gap fixtures add
40 point + 40 SE checks. Total 60 official point and 60 official SE checks.

Limitations: not full certification. No large Monte Carlo; no production-scale
optimization; no full nonabsorbing, PMD, survey-weight, RCS, extra-FE, rank-deficient,
or arbitrary treatment-missingness support. Cluster-by-time nesting VCE degree
conventions outside tested ordinary unit clusters have not been oracle-certified.
Stata process isolation was sequential; concurrent licensed sessions not tested.
