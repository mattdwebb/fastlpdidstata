clear all
set more off
capture log close
log using inference_test_results.log, text replace
adopath ++ "."
use tests/panel.dta, clear
matrix A=(0,3,3,3,3,5,5,5,5,0,0,0 \ 0,0,3,3,3,3,5,5,5,5,0,0)
fastlpdid y, unit(i) time(t) gvar(g) window(-2 2) inference(jackknife ri_beta ri_t_native) assignments(A)
scalar bb=e(estimate)
scalar ns=e(native_se)
scalar js=e(jk_se)
scalar pb=e(ri_beta_p)
scalar pt=e(ri_t_native_p)
assert e(jk_df)==11
matrix dd=e(jk_deletions)
assert rowsof(dd)==12
use tests/r_inference.dta, clear
assert abs(estimate-bb)<1e-10
assert abs(native_se-ns)<1e-10
assert abs(jk_se-js)<1e-10
assert abs(ri_beta-pb)<1e-10
assert abs(ri_t_native-pt)<1e-10
di "FASTLPDID_STATA_INFERENCE_PASS"
log close
exit, clear
