clear all
set more off
capture log close
log using oracle_test_results.log, text replace
sysdir set PLUS "oracle/ado"
sysdir set PERSONAL "oracle/ado"
adopath ++ "."
adopath ++ "oracle"
use tests/panel.dta, clear
xtset i t
gen double lx=L.x
tempname handle
postfile `handle' str12 weighting byte covariates int horizon double estimate double se using tests/stata_oracle.dta, replace
foreach w in variance reweighted {
    foreach cv in 0 1 {
        local rw
        if "`w'"=="reweighted" local rw rw
        local ctl
        local cov
        if `cv' {
            local ctl controls(lx)
            local cov covariates(x)
        }
        lpdid y, unit(i) time(t) treat(d) pre_window(2) post_window(2) pooled nograph `rw' `ctl'
        matrix pp=e(pooled_results)
        matrix hh=e(results)
        matrix list pp
        matrix list hh
        fastlpdid y, unit(i) time(t) treatment(d) window(-2 2) weighting(`w') `cov'
        scalar nb=e(estimate)
        scalar ns=e(native_se)
        matrix nn=e(horizons)
        assert abs(nb-pp[2,1])<1e-9
        assert abs(ns-pp[2,2])<1e-9
        post `handle' ("`w'") (`cv') (999) (pp[2,1]) (pp[2,2])
        forvalues r=1/5 {
            if `r'!=2 {
                assert abs(nn[`r',2]-hh[`r',1])<1e-9
                assert abs(nn[`r',3]-hh[`r',2])<1e-9
                local h=`r'-3
                post `handle' ("`w'") (`cv') (`h') (hh[`r',1]) (hh[`r',2])
            }
        }
    }
}
postclose `handle'
use tests/stata_oracle.dta, clear
format estimate se %24.17g
export delimited using tests/stata_oracle.csv, replace datafmt
di "FASTLPDID_OFFICIAL_ORACLE_PASS"
log close
exit, clear
