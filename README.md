# Local-projections DiD in Stata

Estimate policy effects at chosen horizons using a genuine panel of schools.

The Dube–Girardi–Jordà–Taylor LP-DiD estimator compares outcome changes around policy adoption with changes for eligible controls. It estimates horizon-specific effects and a pooled effect. Reweighted pooling is the default; variance weighting is an explicit alternative.

## Install

Requires **Stata 16 or later**. No R, Python, or additional Stata estimation
packages are required to use this command.

Install the command directly from GitHub:

```stata
net install fastlpdid, from("https://raw.githubusercontent.com/mattdwebb/fastlpdidstata/main/") replace
discard
which fastlpdid
```

For a source-only trial, `adopath ++ "."` makes this checkout available for the
current session. Use one known version and record it for reproducible work.

## A school-policy example

The included fictional data contain 96 school-year observations: 12 schools followed from 2010 through 2017.
Four schools adopt in 2013, four in 2015, and four never adopt. The outcome is
`test_score`, measured in score points. `baseline_score` is a predetermined
covariate. Schools are the independent inference clusters.

The data generator adds a two-point policy effect, alongside school, year,
covariate, and varying disturbance components. This is a software teaching
example, not evidence about an actual policy. The [data guide](examples/DATA_GUIDE.md)
defines every variable and links the reproducible generator.

Run from the extracted repository folder:

```stata
clear
set type double
insheet using "examples/school_panel.csv", comma names clear
fastlpdid test_score, unit(school) time(year) gvar(first_treat) window(0 1) covariates(baseline_score) inference(native jackknife)
display e(estimate)
display e(jk_se)
```

The executed example gives:

| Estimate (score points) | Native clustered SE | Jackknife SE |
|---:|---:|---:|
| 1.601279 | 0.331234 | 0.365699 |

The estimated effect is about **1.60 score points**. The jackknife
standard error is **0.37 points**, calculated by leaving each of the twelve
schools out in turn. It describes uncertainty around this estimate; it does not
change the point estimate or establish the identifying assumptions of a policy
study. Standalone jackknife p-values use Student t with 11 degrees of freedom here.

`do examples/school_policy.do` runs this example, prints a compact table, and
saves `examples/school_results.csv`. Its [recorded output](examples/expected_output.txt)
is checked against the release example. Small last-digit differences may depend
on numerical libraries; they should not change the substantive result.

## Adapt it to your data

Use zero for never-treated schools and repeat each school’s first treatment year on all its rows. Treatment must be absorbing and common to observations in a school-year.

Covariates go in `covariates(baseline_score)`; omit that option for an unadjusted specification. Use numeric covariates whose timing makes sense for your study;
the command does not make post-treatment controls causally appropriate.

Request `inference(native)` for the estimator’s analytic clustered calculation, `inference(jackknife)` for delete-school inference, or `inference(native jackknife)` to retain both. Native and jackknife standard errors are stored separately.

When a deletion removes a cohort, the jackknife recomputes the estimator’s
supported aggregate and its weights using the surviving estimable effects.
All cluster deletions must yield a scalar estimate. A genuinely unidentified
deletion is reported as a failure; it is never silently omitted.

Coefficient randomization inference and the supported studentized variants use
alternative cohort assignments while preserving treatment counts. These require
more input than the introductory example; see [inference options](INFERENCE.md).
Changing the treatment assignment must also change any treatment-dependent fits.

## Supported scope

LP-DiD requires a genuine panel and an explicit event window. The example uses the same schools over time, not pupils linked across repeated cross sections. Covariates are lagged one period by default. Repeated-cross-section equations, observation weights, and nonabsorbing treatment are unsupported.

Different DiD estimators can target different weighted averages when effects
vary. Comparing their coefficients does not by itself show that one is correct.
Unsupported options return explicit errors rather than substituting another model.

## Help, methods, and validation

- [Estimator definition and references](METHOD_SPECIFICATION.md)
- [Inference definitions](INFERENCE.md)
- [Technical development reference](docs/DEVELOPMENT_REFERENCE.md)
- [Runnable example](examples/school_policy.do) and [data definitions](examples/DATA_GUIDE.md)
- [License](LICENSE)

Use the repository’s Issues page for a reproducible problem report. Include your
software version, a small artificial dataset, and the exact command. Do not post
confidential research data. For developers, `tests/run_release_tests.ps1` runs the
native Stata checks and records their status with a licensed Stata installation.
