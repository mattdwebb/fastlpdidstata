{smcl}
{title:fastlpdid (development reference implementation)}
{p}fastlpdid outcome, unit(id) time(t) treatment(d) window(0 3)
{p}Options: gvar(g) instead of treatment(d); cluster(c); covariates(x1 x2);
weighting(reweighted|variance); covlag(1|0). Default weighting is reweighted.
Covariates default to t-1, not treatment-period values.
{p}Only absorbing-treatment genuine panels are supported. Integer calendar time
uses a period of one; gaps are not compressed. Cohort zero/missing means never.
No default event window. Pooled outcome is mean of post-window outcomes minus
the last pre-entry outcome. e(horizons) contains horizon results.
{p}Native Mata point/cluster-sandwich implementation. No Python, R or lpdid calls.
Inference options: inference(native jackknife ri_beta ri_t_native ).
 is an alias for . RI requires assignments(matrix), one row per
unique alternative assignment, columns in sorted cluster order. Zero means never.
All cohort counts must be preserved; the observed assignment must be excluded.
Only null(0) is supported. Stata does not generate assignments in this build.
{p}Validated against official lpdid 1.1.0 on deterministic panel fixtures;
this is not a claim of full certification. See VALIDATION.md.
