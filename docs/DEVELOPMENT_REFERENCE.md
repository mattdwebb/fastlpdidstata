> Historical development notes retained for provenance. For current behavior, installation and examples, use the [main user guide](../README.md). Old audit statements below may describe superseded behavior.

# fastlpdidstata

Audit correction (2026-09-30): `e(native_df)` now exposes the native t(G-1)
reference, or missing for the covariate regression-adjustment normal reference.
Standalone jackknife rejects deletions that change pooled cohort/time support.
Its cluster scope remains all input inference clusters, reported in `e(jk_N_clust)`
and `e(jk_cluster_scope)` separately from native `e(N_clust)`.

Native Stata/Mata reference DGJT LP-DiD. No lpdid, R or Python runtime dependency.
MIT original code; frozen test-only upstream oracle under oracle/ remains GPL-3.
Tested development subset, not full-feature certification.

Install with `net install fastlpdid, from("<checkout-directory>")` or add the
checkout to adopath. The .pkg installs only native code/help, not oracle files.

```stata
adopath ++ "."
use tests/panel.dta, clear
fastlpdid y, unit(i) time(t) gvar(g) window(-2 2) inference(native jackknife)
matrix list e(horizons)
```

Run test.do and inference_test.do. Official validation additionally needs the
test-only dependencies installed by install_oracle.do into oracle/ado (not global).
Run oracle_test.do and oracle_edges.do; require explicit PASS markers in the
distinct *_results.log, not merely a zero OS exit code. Logs are ignored because
Stata startup logs can contain license information.
