*! fastlpdid 0.1.0 29sep2026 -- native reference implementation
program define fastlpdid, eclass sortpreserve
    version 16
    syntax varname(numeric), Unit(varname numeric) Time(varname numeric) ///
        WINDOW(numlist integer min=2 max=2) [TREATment(varname numeric) ///
        GVAR(varname numeric) CLuster(varname) COVariates(varlist numeric) ///
        WEIGHTing(string) COVLag(integer 1) DATAMode(string) ///
        INference(string) ASSIGNments(name) NULL(real 0)]
    if `null'!=0 {
        di as error "UNSUPPORTED: only null(0); nonzero sharp-null adjustment not implemented"
        exit 198
    }
    if "`cluster'"=="" local cluster `unit'
    local opts `varlist', unit(`unit') time(`time') window(`window') cluster(`cluster') covlag(`covlag')
    if "`weighting'"!="" local opts `opts' weighting(`weighting')
    if "`datamode'"!="" local opts `opts' datamode(`datamode')
    if "`covariates'"!="" local opts `opts' covariates(`covariates')
    local base `opts'
    if "`gvar'"!="" local opts `opts' gvar(`gvar')
    if "`treatment'"!="" local opts `opts' treatment(`treatment')
    if "`inference'"=="" local inference native
    local wantjk=0
    local wantri=0
    foreach m of local inference {
        if !inlist("`m'","point","native","jackknife","ri_beta","ri_t_native") error 198
        if inlist("`m'","jackknife") local wantjk=1
        if inlist("`m'","ri_beta","ri_t_native") local wantri=1
    }
    quietly _fastlpdid_point `opts'
    tempname saved obsb obss jse jkvals
    estimates store `saved'
    scalar `obsb'=e(estimate)
    scalar `obss'=e(native_se)
    tempvar cg
    quietly egen long `cg'=group(`cluster')
    quietly summarize `cg', meanonly
    local G=r(max)
    if `wantjk' {
        quietly _fastlpdid_jk, cl(`cg') options(`"`opts'"')
        scalar `jse'=r(se)
        matrix `jkvals'=r(deletions)
    }
    tempname p_beta p_t draws observed
    if `wantri' {
        if "`assignments'"=="" {
            di as error "Stata RI requires supplied assignments(matrix), rows by sorted clusters"
            exit 198
        }
        tempvar cohort assigned
        if "`gvar'"!="" quietly gen double `cohort'=`gvar'
        else {
            quietly egen double `cohort'=min(cond(`treatment'==1,`time',.)), by(`unit')
        }
        quietly replace `cohort'=0 if missing(`cohort')
        quietly bysort `cg': assert `cohort'==`cohort'[1]
        matrix `observed'=J(1,`G',.)
        forvalues g=1/`G' {
            quietly summarize `cohort' if `cg'==`g', meanonly
            matrix `observed'[1,`g']=r(mean)
        }
        mata: _fastlpdid_checkassign("`assignments'","`observed'")
        local R=rowsof(`assignments')
        matrix `draws'=J(`R',2,.)
        local rb=0
        local rt=0
        quietly gen double `assigned'=.
        forvalues r=1/`R' {
            forvalues g=1/`G' {
                quietly replace `assigned'=`assignments'[`r',`g'] if `cg'==`g'
            }
            local ropts `base' gvar(`assigned')
            quietly _fastlpdid_point `ropts'
            local b=e(estimate)-`null'
            local s=e(native_se)
            matrix `draws'[`r',1]=`b'
            if abs(`b')>=abs(`obsb'-`null')-64*c(epsdouble)*max(1,abs(`b'),abs(`obsb'-`null')) local ++rb
            if strpos("`inference'","ri_t_native") {
                if missing(`s') | `s'<=c(epsdouble) | missing(`obss') | `obss'<=c(epsdouble) error 498
                matrix `draws'[`r',2]=`b'/`s'
                if abs(`b'/`s')>=abs((`obsb'-`null')/`obss')-64*c(epsdouble)*max(1,abs(`b'/`s'),abs((`obsb'-`null')/`obss')) local ++rt
            }
        }
        scalar `p_beta'=(1+`rb')/(`R'+1)
        scalar `p_t'=(1+`rt')/(`R'+1)
    }
    estimates restore `saved'
    estimates drop `saved'
    if `wantjk' {
        ereturn scalar jk_se=`jse'
        ereturn scalar jk_N_clust=`G'
        ereturn local jk_cluster_scope "all input inference clusters"
        ereturn local jk_support "recomputed pooled cohort/time support"
        ereturn scalar jk_df=`G'-1
        ereturn scalar jk_p=2*ttail(`G'-1,abs((`obsb'-`null')/`jse'))
        ereturn matrix jk_deletions=`jkvals'
    }
    if `wantri' {
        ereturn matrix ri_draws=`draws'
        if strpos("`inference'","ri_beta") ereturn scalar ri_beta_p=`p_beta'
        if strpos("`inference'","ri_t_native") ereturn scalar ri_t_native_p=`p_t'
    }
    ereturn local cmd "fastlpdid"
    ereturn display
end

program define _fastlpdid_jk, rclass
    version 16
    syntax, CL(varname) OPTIONS(string)
    quietly _fastlpdid_point `options'
    tempname fullb
    scalar `fullb'=e(estimate)
    quietly levelsof `cl', local(groups)
    local G : word count `groups'
    tempname vals
    matrix `vals'=J(`G',1,.)
    local j=0
    foreach g of local groups {
        local ++j
        preserve
        quietly keep if `cl'!=`g'
        capture quietly _fastlpdid_point `options'
        local rc=_rc
        if !`rc' matrix `vals'[`j',1]=e(estimate)
        restore
        if `rc' {
            di as error "Jackknife deletion `g' failed; all deletions are required"
            exit `rc'
        }
    }
    tempname v
    mata: st_numscalar("`v'", (`G'-1)/`G'*sum((st_matrix("`vals'"):-mean(st_matrix("`vals'"))):^2))
    tempname scale
    mata: st_numscalar("`scale'",max((1,abs(st_numscalar("`fullb'")),max(abs(st_matrix("`vals'"))))))
    if missing(`v') | sqrt(`v')<=64*2.220446049250313e-16*`scale' exit 498
    return scalar se=sqrt(`v')
    return matrix deletions=`vals'
end

program define _fastlpdid_point, eclass sortpreserve
    version 16
    syntax varname(numeric), Unit(varname numeric) Time(varname numeric) ///
        WINDOW(numlist integer min=2 max=2) [TREATment(varname numeric) ///
        GVAR(varname numeric) CLuster(varname) COVariates(varlist numeric) ///
        WEIGHTing(string) COVLag(integer 1) DATAMode(string)]
    if "`weighting'" == "" local weighting reweighted
    if !inlist("`weighting'", "reweighted", "variance") error 198
    if !inlist("`datamode'", "", "auto", "panel") {
        di as error "UNSUPPORTED: genuine panel required"
        exit 198
    }
    if !inlist(`covlag',0,1) error 198
    if "`treatment'`gvar'" == "" error 198
    if "`cluster'" == "" local cluster `unit'
    local low : word 1 of `window'
    local high : word 2 of `window'
    if `low'>`high' | `high'<0 error 198
    local start = max(0,`low')
    preserve
    quietly {
        isid `unit' `time'
        assert !missing(`unit',`time')
        assert `time' == floor(`time')
        assert !missing(`cluster')
        bysort `unit' (`time'): assert `cluster' == `cluster'[1]
        xtset `unit' `time', delta(1)
        tempvar d g dt c yy ok w nt nc nn
        egen long `c' = group(`cluster')
        if "`gvar'" != "" {
            gen double `g' = `gvar'
            replace `g' = . if `g'==0
            bysort `unit' (`time'): assert `g'==`g'[1]
            assert `g'==floor(`g') if !missing(`g')
            gen byte `d' = `time'>=`g' & !missing(`g')
            if "`treatment'" != "" assert `d' == `treatment'
        }
        else {
            gen byte `d' = `treatment'
            assert inlist(`d',0,1)
            bysort `unit' (`time'): assert `d'>=`d'[_n-1] if _n>1
            by `unit': assert `time'==`time'[_n-1]+1 if _n>1 & `d'==1 & `d'[_n-1]==0
        }
        sort `unit' `time'
        gen double `dt' = D.`d'
        local xs
        foreach v of local covariates {
            tempvar x
            if `covlag'==1 gen double `x' = L.`v'
            else gen double `x' = `v'
            local xs `xs' `x'
        }
        local xopt
        if "`xs'"!="" local xopt xs(`xs')
        tempname results pooled ans b V
        local rows = `high'-`low'+1
        matrix `results' = J(`rows',5,.)
        forvalues h = `low'/`high' {
            local row = `h'-`low'+1
            matrix `results'[`row',1] = `h'
            if `h'==-1 {
                matrix `results'[`row',2] = 0
                continue
            }
            if `h'<0 {
                local ah = -`h'
                gen double `yy' = L`ah'.`varlist'-L.`varlist'
            }
            else if `h'==0 gen double `yy' = `varlist'-L.`varlist'
            else gen double `yy' = F`h'.`varlist'-L.`varlist'
            local end = max(0,`h')
            if `end'==0 gen byte `ok' = `dt'==1 | `d'==0
            else gen byte `ok' = `dt'==1 | F`end'.`d'==0
            markout `ok' `yy' `dt' `xs' `c'
            capture _fastlpdid_reg `yy', d(`dt') time(`time') cl(`c') ok(`ok') `xopt' weighting(`weighting')
            if !_rc {
                matrix `ans' = r(result)
                forvalues j=1/4 {
                    matrix `results'[`row',`j'+1] = `ans'[1,`j']
                }
            }
            drop `yy' `ok'
        }
        gen double `yy' = 0
        forvalues h=`start'/`high' {
            if `h'==0 replace `yy' = `yy'+`varlist'
            else replace `yy' = `yy'+F`h'.`varlist'
        }
        replace `yy' = `yy'/(`high'-`start'+1)-L.`varlist'
        if `high'==0 gen byte `ok' = `dt'==1 | `d'==0
        else gen byte `ok' = `dt'==1 | F`high'.`d'==0
        markout `ok' `yy' `dt' `xs' `c'
        _fastlpdid_reg `yy', d(`dt') time(`time') cl(`c') ok(`ok') `xopt' weighting(`weighting')
        matrix `pooled' = r(result)
        matrix `b' = `pooled'[1,1]
        matrix colnames `b' = pooled
        matrix `V' = (`pooled'[1,2]^2)
        matrix colnames `V' = pooled
        matrix rownames `V' = pooled
        matrix colnames `results' = horizon estimate se N clusters
        local N = `pooled'[1,3]
        local G = `pooled'[1,4]
        levelsof `time' if `ok', local(pooledtimes)
        levelsof `time' if `ok' & `dt'==1, local(pooledcohorts)
    }
    restore
    ereturn post `b' `V', obs(`N')
    ereturn scalar estimate = `pooled'[1,1]
    ereturn scalar native_se = `pooled'[1,2]
    ereturn scalar N_clust = `G'
    ereturn local pooled_times "`pooledtimes'"
    ereturn local pooled_cohorts "`pooledcohorts'"
    ereturn scalar native_df=.
    if "`weighting'"=="variance" | "`covariates'"=="" {
        ereturn scalar df_r=`G'-1
        ereturn scalar native_df=`G'-1
        ereturn scalar native_p=2*ttail(`G'-1,abs(`pooled'[1,1]/`pooled'[1,2]))
    }
    else ereturn scalar native_p=2*normal(-abs(`pooled'[1,1]/`pooled'[1,2]))
    ereturn matrix horizons = `results'
    ereturn local weighting "`weighting'"
    ereturn local target_label "`weighting' pooled LP-DiD ATT [`start',`high']"
    ereturn local validation "Official lpdid 1.1.0 deterministic fixture; not full certification"
    ereturn local cmd "fastlpdid"
    ereturn display
end

program define _fastlpdid_reg, rclass sortpreserve
    version 16
    syntax varname, D(varname) Time(varname) CL(varname) OK(varname) ///
        WEIGHTing(string) [XS(varlist)]
    tempvar n nt w
    quietly {
        bysort `time': egen double `n' = total(`ok')
        by `time': egen double `nt' = total(`ok'*`d')
        gen double `w' = 1
        if "`weighting'"=="reweighted" & "`xs'"!="" {
            replace `ok' = 0 if `nt'==`n'
        }
        else if "`weighting'"=="reweighted" {
            replace `ok' = 0 if `nt'==0 | `nt'==`n'
            replace `w' = `n'/(`n'-`nt') if `ok'
        }
        else replace `ok' = 0 if `n'<=1
    }
    tempname ans
    local ra = "`weighting'"=="reweighted" & "`xs'"!=""
    mata: _fastlpdid_solve("`varlist'","`d'","`time'","`cl'","`xs'","`w'","`ok'",`ra',"`ans'")
    if missing(`ans'[1,1]) {
        di as error "IDENTIFICATION_FAILURE: singular design or no clean comparison"
        exit 498
    }
    return matrix result = `ans'
end

mata:
void _fastlpdid_checkassign(string scalar an, string scalar on)
{
    real matrix A
    real rowvector O
    real scalar i,j
    A=st_matrix(an); O=st_matrix(on)
    if(cols(A)!=cols(O) | rows(A)<1 | hasmissing(A)) _error(198)
    for(i=1;i<=rows(A);i++) {
        if(any(sort(A[i,]',1):!=sort(O',1)) | all(A[i,]:==O)) _error(198)
        for(j=1;j<i;j++) if(all(A[i,]:==A[j,])) _error(198)
    }
}

void _fastlpdid_solve(string scalar yn, string scalar dn, string scalar tn,
    string scalar cn, string scalar xn, string scalar wn, string scalar sn,
    real scalar ra, string scalar result)
{
    real matrix X, FE, B, X0
    real colvector y,d,t,c,w,times,clusters,b,r,s,cs,ix
    real rowvector m
    real scalar j,N,G,K,est,se,correction
    y=st_data(.,yn,sn); d=st_data(.,dn,sn); t=st_data(.,tn,sn)
    c=st_data(.,cn,sn); w=st_data(.,wn,sn); N=rows(y)
    st_matrix(result,J(1,4,.))
    if(N==0 | sum(d)==0 | sum(1:-d)==0) return
    times=uniqrows(sort(t,1)); clusters=uniqrows(sort(c,1)); G=rows(clusters)
    FE=J(N,rows(times),0)
    for(j=1;j<=rows(times);j++) FE[,j]=(t:==times[j])
    X=FE
    if(strlen(xn)) X=X,st_data(.,tokens(xn),sn)
    if(ra) {
        ix=selectindex(d:==0); X0=X[ix,]
        if(rank(X0)<cols(X0)) return
        B=invsym(cross(X0,X0)); b=B*cross(X0,y[ix])
        r=y-X*b; ix=selectindex(d:==1); est=mean(r[ix]); m=mean(X[ix,])
        s=(d:*(r:-est)-(1:-d):*r:*(X*B*m')*sum(d))/sum(d)
        correction=1
    }
    else {
        X=d,X; K=cols(X)
        if(rank(X)<K) return
        B=invsym(cross(X,w,X)); b=B*cross(X,w,y); est=b[1]
        s=(X*B[,1]):*w:*(y-X*b)
        correction=G/(G-1)*(N-1)/(N-K)
    }
    cs=J(G,1,0)
    for(j=1;j<=G;j++) cs[j]=sum(select(s,c:==clusters[j]))
    se=sqrt(correction*cross(cs,cs))
    st_matrix(result,(est,se,N,G))
}
end
