clear all
set more off
capture log close
log using oracle_edges_results.log, text replace
sysdir set PLUS "oracle/ado"
sysdir set PERSONAL "oracle/ado"
adopath ++ "."
adopath ++ "oracle"
forvalues mode=1/2 {
    use tests/panel.dta, clear
    if `mode'==1 replace y=. if i==1 & t==4
    if `mode'==2 drop if i==1 & t==2
    xtset i t
    gen double lx=L.x
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
            di "EDGE CASE mode=`mode' weight=`w' covariates=`cv'"
            lpdid y, unit(i) time(t) treat(d) pre_window(2) post_window(2) pooled nograph `rw' `ctl'
            matrix pp=e(pooled_results)
            matrix hh=e(results)
            quietly fastlpdid y, unit(i) time(t) gvar(g) window(-2 2) weighting(`w') `cov'
            di "Pooled differences: " %24.17g e(estimate)-pp[2,1] " SE: " %24.17g e(native_se)-pp[2,2]
            assert abs(e(estimate)-pp[2,1])<1e-9
            assert abs(e(native_se)-pp[2,2])<1e-9
            matrix nn=e(horizons)
            forvalues r=1/5 {
                if `r'!=2 {
                    assert abs(nn[`r',2]-hh[`r',1])<1e-9
                    assert abs(nn[`r',3]-hh[`r',2])<1e-9
                }
            }
        }
    }
}
di "FASTLPDID_OFFICIAL_EDGES_PASS"
log close
exit, clear
