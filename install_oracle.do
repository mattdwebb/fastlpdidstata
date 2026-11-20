clear all
set more off
capture log close
log using install_oracle_results.log, text replace
capture mkdir oracle/ado
sysdir set PLUS "oracle/ado"
sysdir set PERSONAL "oracle/ado"
ssc install ftools
ssc install require
ssc install reghdfe
ssc install boottest
ssc install listreg
di "LPDID_ORACLE_DEPS_READY"
log close
exit, clear
