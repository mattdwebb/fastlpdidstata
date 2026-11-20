clear all
set more off
capture set processors 1
capture log close
log using test_results.log, text replace
adopath ++ "."
use tests/panel.dta, clear
foreach w in variance reweighted {
    foreach cv in 0 1 {
        local xopt
        if `cv' local xopt covariates(x)
        fastlpdid y, unit(i) time(t) gvar(g) window(-2 2) weighting(`w') `xopt'
        matrix hh=e(horizons)
        scalar bp=e(estimate)
        scalar sp=e(native_se)
        preserve
        use tests/r_port.dta, clear
        keep if weighting=="`w'" & covariates==`cv'
        assert abs(estimate-bp)<1e-9 if horizon==999
        assert abs(se-sp)<1e-9 if horizon==999
        foreach h in -2 0 1 2 {
            local r=`h'+3
            assert abs(estimate-hh[`r',2])<1e-9 if horizon==`h'
            assert abs(se-hh[`r',3])<1e-9 if horizon==`h'
        }
        restore
    }
}
di "FASTLPDID_STATA_CROSS_LANGUAGE_PASS"
log close
exit, clear
