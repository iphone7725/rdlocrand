
/*******************************************************************************

* Window selection for randomization inference in RD

*!version 0.3 13-Mar-2018

Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare

*******************************************************************************/

version 13

capture program drop rdwinselect
program define rdwinselect, rclass sortpreserve

	syntax varlist (min=1 numeric) [if] [in] [, Cutoff(real 0) ///
												obsmin(numlist max=1) ///												
												wmin(numlist max=1) ///
												wobs(numlist max=1) ///
												wstep(numlist max=1) ///
												NWindows(real 10) ///
												STATistic(string) ///
												APPROXimate ///
												p(integer 0) ///
												evalat(string) ///
												kernel(string) ///
												reps(integer 1000) ///
												seed(integer 666) ///
												level(real .15) ///
												plot ///
												graph_options(string) ///
												genwvars(string) ///
												obsstep(numlist max=1) ]
	
	tokenize `varlist'
	local runv_aux "`1'"
	mac shift 1
	local covariates "`*'"
	marksample touse, novarlist
	marksample touse1
	
	tempvar treated runvar
	qui gen double `runvar' = `runv_aux' - `cutoff' if `touse'
	qui gen `treated' = `runvar' >= 0 if `touse'
	sort `runvar'
	
	quietly summarize `runv_aux' if `touse'
	if(`cutoff' <= r(min) | `cutoff' >= r(max)) {
		display as error "cutoff must be within the range of running variable"
		exit 125
	}
	
	if `p'>0 & "`approximate'"!=""&"`statistic'"!=""&"`statistic'"!="ttest"&"`statistic'"!="diffmeans"{
		di as error "approximate and p>1 can only be combined with diffmeans"
		exit 198
	}
	
	if "`evalat'"!=""&"`evalat'"!="means"&"`evalat'"!="cutoff"{
		di as error "evalat only admits means or cutoff"
		exit 198
	}
	
	if "`kernel'"!=""&"`kernel'"!="uniform"&"`kernel'"!="triangular"&"`kernel'"!="epan"{
		di as error "`kernel' not a valid kernel"
		exit 198
	}
	if "`kernel'"!="" & "`evalat'"!="" &"`evalat'"!="cutoff"{
		di as error "kernel can only be combined with evalat(cutoff)"
		exit 198
	}
	if "`statistic'"!="ttest"&"`statistic'"!=""|"`statistic'"!="diffmeans"{
		if "`kernel'"!=""&"`kernel'"!="uniform"{
			di as error "kernel only allowed for diffmeans"
			exit 198
		}
	}
	
	if `seed'==666{	
		set seed 666
	}
	else if `seed'>=0{
		set seed `seed'
	}
	else if `seed'!=-1{
		di as error "seed must be a positive integer (or -1 for system seed)"
		exit 198
	}
	else {
		local seed c(seed)
	}
	
	if "`genwvars'"!=""{
		tokenize "`genwvars'", parse(",")
		local vars `1'
		local rep_op `3'
		tokenize `vars'
		local wid_var `1'
		local wlength_var `2'
		if "`rep_op'"=="replace"{
			capture drop `wid_var'
			capture drop `wlength_var'
		}
		if "`rep_op'"!="replace" & "`rep_op'"!=""{
			di as error "unrecognized `rep_op'"
			exit 198
		}
	}
	
	tempvar Wid Wlength dups
	qui bysort `runvar': gen `dups' = _N
	qui gen double `Wid' = .
	qui gen double `Wlength' = .


********************************************************************************
** Define initial window

/* Ensures that the smallest window contains at least obsmin observations
on each side of the cutoff (default obsmin=10) */
********************************************************************************
	
	if "`wmin'"!=""{
		if "`obsmin'"!=""{
			di as error "Cannot specify both obsmin and wmin"
			exit 184
		}
		else{
			local obsmin = 10
		}
	}
	else {
		qui count if `treated'==0
		local posl = r(N)
		local posr = r(N) + 1

		if "`obsmin'"!=""{
			qui mata: findwobs(`obsmin',1,`posl',`posr',"`runvar'","`dups'","`Wlength'")
			local wmin = wlength
		}
		else {
			local obsmin = 10
			mata: findwobs(10,1,`posl',`posr',"`runvar'","`dups'","`Wlength'")
			local wmin = wlength
		}
	}
	
	qui replace `Wid' = .
	qui replace `Wlength' = .

	
********************************************************************************
** Define step

/* Default is obsstep=2 for backward compatibility. This calculates a fixed 
length increment to add at least obsstep observations at each side of the cutoff
for all windows. The use of wstep and wobs is recommended. */
********************************************************************************

	if "`wstep'"=="" & "`wobs'"==""{
		if "`obsstep'"==""{
			local obsstep = 2
		}
		tempvar control 
		qui gen `control' = 1 - `treated'
		mata: findstep(`obsmin',`obsstep',10,"`runvar'","`treated'","`control'")
		local wstep = step
		local wmax = `wmin' + `wstep'*(`nwindows'-1)
		local wlist "`wmin'(`wstep')`wmax'"
	}
	
	else if "`wstep'"!=""{
		if "`wobs'"!=""{
			di as error "Cannot specify both wobs and wstep"
			exit 184
		}
		if "`obsstep'"!=""{
			di as error "Cannot specify both obsstep and wstep"
			exit 184
		}
		local wmax = `wmin' + `wstep'*(`nwindows'-1)
		local wlist "`wmin'(`wstep')`wmax'"
	}
	
	else if "`wobs'"!=""{
		qui count if `treated'==0 & `touse'
		local pos0 = r(N)		
		qui count if `runvar'>=-`wmin' & `runvar'<0 & `touse'
		local posl = `pos0' - r(N)
		qui count if `runvar'>=0 & `runvar'<=`wmin' & `touse'
		local posr = `pos0' + 1 + r(N)
		mata: findwobs(`wobs',`nwindows'-1,`posl',`posr',"`runvar'","`dups'","`Wlength'")
		qui levelsof `Wlength', local(wlist)
		local wlist "`wmin' `wlist'"
	}

	
	
********************************************************************************
** Sample sizes
********************************************************************************
	
	tempvar treat_aux contr_aux
	qui{
		gen `treat_aux' = `runvar' if `treated'==1
		count if `treated'==1 & `touse'
		local nt=r(N)
		_pctile `treat_aux' if `touse', p(1 5 10 20)
		local p1t=round(r(r1),.00001)
		local p5t=round(r(r2),.00001)
		local p10t=round(r(r3),.00001)
		local p20t=round(r(r4),.00001)
		count if `treated'==1 & `runvar'<=`p1t' & `touse'
		local nt1=r(N)
		count if `treated'==1 & `runvar'<=`p5t' & `touse'
		local nt5=r(N)
		count if `treated'==1 & `runvar'<=`p10t' & `touse'
		local nt10=r(N)
		count if `treated'==1 & `runvar'<=`p20t' & `touse'
		local nt20=r(N)
		
		gen `contr_aux'=abs(`runvar') if `treated'==0
		count  if `treated'==0 & `touse'
		local nc=r(N)
		_pctile `contr_aux' if `touse', p(1 5 10 20)
		local p1c=r(r1)
		local p5c=r(r2)
		local p10c=r(r3)
		local p20c=r(r4)
		count if `treated'==0 & `runvar'>=-`p1c' & `touse'
		local nc1=r(N)
		count if `treated'==0 & `runvar'>=-`p5c' & `touse'
		local nc5=r(N)
		count if `treated'==0 & `runvar'>=-`p10c' & `touse'
		local nc10=r(N)
		count if `treated'==0 & `runvar'>=-`p20c' & `touse'
		local nc20=r(N)		
		
	}

	
********************************************************************************	
** Define statistic
********************************************************************************

	if "`statistic'"=="" | "`statistic'"=="diffmeans" | "`statistic'"=="ttest"{
		local command "diffmeans"
		local stat_rdri "stat(ttest)"
	}
	else {
		if "`statistic'"=="ksmirnov"{
			local command "ksmirnov"
			local stat_rdri "stat(ksmirnov)"
		}
		else if "`statistic'"=="ranksum"{
			local command "ranksum"
			local stat_rdri "stat(ranksum)"
		}
		else if "`statistic'"=="hotelling"{
			local command "hotelling"
			if `p'>0{
				di as error "p()>0 not allowed for Hotelling statistic"
				exit 198
			}
		}
		else{
			di as error "`statistic' not a valid statistic"
			exit 198
		}
	}

	if "`approximate'"==""{
		local di_method "rdrandinf"
		local di_reps=`reps'
	}
	else {
		local di_method "approximate"
		local di_reps " ."
	}
	
	if `p'==0{
		local model "."
	}
	else if `p'==1{
		local model "linear"
	}
	else {
		local model "polynomial"
	}
	
	if "`kernel'"=="uniform"|"`kernel'"==""{
		local kernel_disp "uniform"
	}
	
	if "`kernel'"=="triangular" {
		local kernel_disp "triangular"
	}

	if "`kernel'"=="epan" {
		local kernel_disp "Epanechnikov"
	}
	
	
********************************************************************************
** Display upper panel
********************************************************************************

	di _newline
	di as text "Window selection for RD under local randomization"
	di _newline
	disp as text "Cutoff c = " as res %4.2f `cutoff'  	as text	_col(18) " {c |} " _col(19) 	 "Left of c"	_col(33) in gr "Right of c" 	_col(51) 			"Number of obs  = " as res %13.0f `nt'+`nc'
	disp as text "{hline 18}{c +}{hline 23}"                                                                                              		_col(51) 			"Order of poly  = " as res %13.0f `p'
	disp as text "{ralign 17:Number of obs}"    				_col(18) " {c |} " _col(17) as res %9.0f `nc'     _col(34) %9.0f  `nt'  	 	_col(51) as text 	"Kernel type    = "	as res "{ralign 13: `kernel_disp'}"
	disp as text "{ralign 17:1st percentile}"   				_col(18) " {c |} " _col(17) as res %9.0f `nc1'    _col(34) %9.0f  `nt1' 	 	_col(51) as text	"Reps           = " as res %13.0f `di_reps'
	disp as text "{ralign 17:5th percentile}"   				_col(18) " {c |} " _col(17) as res %9.0f `nc5'    _col(34) %9.0f  `nt5' 	 	_col(51) as text	"Testing method = " as res "{ralign 13: `di_method'}"
	disp as text "{ralign 17:10th percentile}"  				_col(18) " {c |} " _col(17) as res %9.0f `nc10'   _col(34) %9.0f  `nt10' 		_col(51) as text	"Balance test   = " as res "{ralign 13: `command'}"
	disp as text "{ralign 17:20th percentile}"  				_col(18) " {c |} " _col(17) as res %9.0f `nc20'   _col(34) %9.0f  `nt20' 


********************************************************************************
** Main results
********************************************************************************

	local matrows ""
	mat Results = J(`nwindows',5,.)

	di _newline
	di as text 			_col(18) " {c |}" 			_col(23) "Bal. test " 	_col(41) "Var. name"				_col(54) "Bin. test "
	di as text _col(1) " Window length /2 {c |}" 	_col(24) "p-value" 		_col(39) "(min p-value)" 		_col(55) "p-value" 	_col(66)" Obs<c " _col(74) " Obs>=c"
	di as text "{hline 18}{c +}{hline 61}"

	local j = 1
	foreach w of numlist `wlist'{

		local wupper = `w'
		local wlower = -`w'
	
		local wuppers: di %4.3f `cutoff' + `wupper'

		local inwindow "`runvar' >= float(`wlower') & `runvar' <= float(`wupper')"
		
		qui replace `Wid' = `j' if `inwindow' & `Wid'==. & `touse'
		qui replace `Wlength' = `w' if `inwindow' & `Wlength'==. & `touse'

		preserve
		qui keep if `inwindow' & `touse1'
		
		* Covariate balance test *
		
		if "`covariates'"!=""{
			local ncov: word count `covariates'

			* Model adjustment 
			
			tempvar Y_adj Y_adj_null kweights
			
			if "`kernel'"=="triangular" {
				local bwt = `w'
				local bwc = `wupper'
				qui gen `kweights' = 1-abs((`cutoff'-`runv_aux')/`bwt') if abs((`cutoff'-`runv_aux')/`bwt')<1 & `treated'==1
				qui replace `kweights' = 1-abs((`cutoff'-`runv_aux')/`bwc') if abs((`cutoff'-`runv_aux')/`bwc')<1 & `treated'==0
				local kweights_opt "[aw = `kweights']"
				local kwrd_opt "weights(`kweights')"
			}

			if "`kernel'"=="epan" {
				local bwt = `wupper'
				local bwc = `wupper'
				qui gen `kweights' = .75*(1-((`cutoff'-`runv_aux')/`bwt')^2) if abs((`cutoff'-`runv_aux')/`bwt')<1 & `treated'==1
				qui replace `kweights' = .75*(1-((`cutoff'-`runv_aux')/`bwc')^2) if abs((`cutoff'-`runv_aux')/`bwc')<1 & `treated'==0
				local kweights_opt "[aw = `kweights']"
				local kwrd_opt "weights(`kweights')"
			}
			
			if `p'>0{
				if "`evalat'"==""|"`evalat'"=="cutoff"{
					local evalr = `cutoff'
					local evall = `cutoff'
				}
				else {
					qui sum `runv_aux' if `treated'==1
					local evalr = r(mean)
					qui sum `runv_aux' if `treated'==0
					local evall = r(mean)
				}
				tempvar r_t r_c resid_l resid_r
				qui gen double `r_t' = `runv_aux'-`evalr'
				qui gen double `r_c' = `runv_aux'-`evall'
				
				foreach cov of varlist `covariates'{
					qui{
					
						forvalues k=1/`p'{
							gen _rpt_`cov'`k'=`r_t'^`k'
						}
						reg `cov' _rpt_`cov'* if `treated'==1 `kweights_opt'
						predict `resid_r' if e(sample), residuals
						gen double _adj_`cov' = `resid_r' + _b[_cons] if e(sample)
						
						forvalues k=1/`p'{
							gen _rpc_`cov'`k'=`r_c'^`k'
						}
						reg `cov' _rpc_`cov'* if `treated'==0 `kweights_opt'
						predict `resid_l' if e(sample), residuals
						replace _adj_`cov' = `resid_l' + _b[_cons] if e(sample)
						
						drop `resid_l' `resid_r'
					}
				}
			}

			* Balance test
			
			if "`statistic'"=="hotelling"{
				if "`approximate'"!=""{	
					qui hotelling `covariates', by(`treated')
					local df1 = r(k)
					local df2 = r(N)-1-r(k)
					local pval_h = 1-F(`df1',`df2',r(T2))
				}
				else{
					qui permute `treated' stat=r(T2), reps(`reps') nowarn nodots: /*
						*/ hotelling `covariates', by(`treated')
					mat Pvals = r(p)'
					local pval_h = Pvals[1,1]
				}
				local xminp "-"
			}
				
			else{
				matrix Pvals=J(`ncov',1,.)
				local row=1	
				
				if "`approximate'"!=""{	
					foreach cov of varlist `covariates'{
						local cov_`row' `cov'
						local dcov_`row' `cov'
							if "`statistic'"=="ttest"|"`statistic'"=="diffmeans"|"`statistic'"==""{
								if `p'==0{
									qui reg `cov' `treated' `kweights_opt', vce(hc2)
									local asy_p = 2*normal(-abs(_b[`treated']/_se[`treated']))
								}
								else {
									if "`evalat'"==""|"`evalat'"=="cutoff"{
										forvalues k=1/`p'{
											gen _rp_`cov'`k'=`runvar'^`k'
										}
										qui reg `cov' `treated'##c.(_rp_`cov'*) `kweights_opt', vce(hc2)	
										local asy_p = 2*normal(-abs(_b[1.`treated']/_se[1.`treated']))
									}
									else {
										qui reg `cov' _runpoly_t_* if `treated'==1
										local a_t = _b[_cons]
										local se_t = _se[_cons]
										qui reg `cov' _runpoly_c_* if `treated'==0
										local a_c = _b[_cons]
										local se_c = _se[_cons]
										local asy_p = 2*normal(-abs(`obs_stat'/sqrt(`se_t'^2+`se_c'^2)))
									}
								}
							}
							else {
								rdrandinf_model `cov' `treated', stat(`statistic') asy
								local asy_p = r(asy_pval)
							}

						mat Pvals[`row',1]=`asy_p'
						local ++row
					}
				}
			
				else {
					local stat_list ""
					if `p'==0{
						foreach cov in `covariates'{
							local cov_`row' `cov'
							local dcov_`row' `cov'
							local stat_list "`stat_list' stat_`row'=r(stat_`row')"
							local ++row
						}
						qui{
							permute `treated' `stat_list', reps(`reps') nowarn nodots: /*
								*/ rdwinselect_allcovs `covariates', treat(`treated') runvar(`runvar') /*
								*/ `stat_rdri' `kwrd_opt'
							mat Pvals = r(p)'
						}
					}
					else {
						foreach cov in `covariates' {
							local cov_`row' _adj_`cov'
							local dcov_`row' `cov'
							local stat_list "`stat_list' stat_`row'=r(stat_`row')"
							local ++row
						}
						qui{
							permute `treated' `stat_list', reps(`reps') nowarn nodots: /*
								*/ rdwinselect_allcovs _adj_*, treat(`treated') runvar(`runvar') /*
								*/ `stat_rdri' `kwrd_opt'
							mat Pvals = r(p)'
						}
					}
				}
				
				mata: minindex(st_matrix("Pvals"),1,minindpi=.,minindpj=.)
				mata: st_numscalar("minp",min(st_matrix("Pvals")))
				mata: st_numscalar("minindp",minindpi[1,1])
				local ind = minindp
				local xminp "`dcov_`ind''"
			}
		
		}
		
		else {
			matrix Pvals = J(1,1,.)
			local xminp "-"
		}

		* Binomial test *
		qui bitest `treated' = 1/2
		mat Results[`j',2] = r(p)
		local pbin = r(p)

		* Number of treated and controls *
		qui count if `treated'==1
		mat Results[`j',4] = r(N)
		local nt = r(N)
		qui count if `treated'==0
		mat Results[`j',3] = r(N)
		local nc = r(N)
	
		* Window length *
		mat Results[`j',5] = `wuppers' - `cutoff'
		if "`covariates'"!=""{
			if "`statistic'"!="hotelling"{
				mat Results[`j',1] = minp
				output_line `wuppers'-`cutoff' minp `xminp' `pbin' `nt' `nc'

			}
			else {
				mat Results[`j',1] = `pval_h'
				output_line `wuppers'-`cutoff' `pval_h' `xminp' `pbin' `nt' `nc'

			}
		}
		else {
			mat Results[`j',1] = .
			output_line `wuppers'-`cutoff' . `xminp' `pbin' `nt' `nc'
		}
		local matrows " `matrows' "`j'""
		
		restore
		local ++j
	}
	
	mat colnames Results = "Bal test" "Bin test" "Obs<c" "Obs>=c" "Length/2"
	mat rownames Results = `matrows'
	
	di _newline as text "Variable used in binomial test (running variable): " as res "`runv_aux'"
	if "`covariates'"!=""{
		local ls_prev = c(linesize)
		set linesize 80
		di as text "Covariates used in balance test: " as res "`covariates'"
		set linesize `ls_prev'
		if "`approximate'"!="" & `p'>0 & ("`statistic'"=="ranksum" | "`statistic'"=="ksmirnov"){
			di as error "warning: asymptotic p-values do not account for outcome model adjustment."
		}
	}
	else {
		di _newline as text "Note: no covariates specified."
	}


********************************************************************************
** Recommended window: largest window before first crossing
********************************************************************************

	if "`covariates'"!=""{
		
		mata: Q=st_matrix("Results")
		mata: reclength(Q[,1],Q[,5],`level')
		
		if reclength==0{
			di _newline as error "Smallest window doesn't pass covariate test. " /*
				*/ "Consider modifying smallest window or reducing the level."
			
			local Nr = .
			local Nl = .
			local Nw = .
			local minp = .
		}
		else {

			di _newline as text "Recommended window is" as res " [" round(`cutoff'-reclength,.001) "; " round(`cutoff'+reclength,.001) "]" /*
				 */ as text " with " as res Results[index,3]+Results[index,4] as text " observations (" as res Results[index,3] /*
				 */ as text " below, " as res Results[index,4] as text " above)."
			local Nw = Results[index,3]+Results[index,4]
			local Nl = Results[index,3]
			local Nr = Results[index,4]
			local minp = Results[index,1]
		}
	}
	else {
		di "Need to specify covariates to find recommended length."
	}	
	
	
********************************************************************************
** Generate plot (if specified)
********************************************************************************

	if "`plot'"!=""{
		if "`covariates'"==""{
			di as error "Cannot draw plot without covariates."
			exit 498
		}
		else {
			preserve
			capture drop _plotvars*
			qui svmat Results, names(_plotvars)
			local xlabels ""
			forvalues i=1(4)`nwindows'{
				local xlabel_`i': di %4.2f _plotvars5[`i']
				local xlabels "`xlabels' `xlabel_`i''"
			}
			if "`graph_options'"==""{
				twoway scatter _plotvars1 _plotvars5, title("Minimum p-value from covariate test") /*
					*/ xtitle("window length / 2") ytitle(P-value) /*
					*/ xlabel(`xlabels') yline(`level', lpattern(shortdash) lcolor(black)) /*
					*/  ysc(r(0)) /*
					*/ note(The dotted line corresponds to p-value=`level')
			}
			else {
				twoway scatter _plotvars1 _plotvars5, `graph_options'
			}
			drop _plotvars*
			restore
		}
	}
	
	
********************************************************************************
** Generate window id variables (if specified)
********************************************************************************

	if "`genwvars'"!=""{
		qui gen double `wid_var' = `Wid'
		qui gen double `wlength_var' = `Wlength'
	}
	
	
********************************************************************************
** Return values
********************************************************************************
	
	return local seed = `seed'
	return matrix results = Results
	
	if "`covariates'"!=""{
		return scalar rec_length = reclength
		return scalar N_right = `Nr'
		return scalar N_left = `Nl'
		return scalar N = `Nw'
		return scalar minp = `minp'
	}
	else {
		return scalar rec_length = .
	}
	
end




********************************************************************************
********************************************************************************
** Auxiliary functions 
********************************************************************************
********************************************************************************


********************************************************************************
** Output display
********************************************************************************

capture program drop output_line
program output_line
	args window cov vname bin nt nc
	display as res %17.3f `window' 	_col(18) as text " {c |}" as result _col(26) %4.3f `cov' " " _col(38) %10s abbrev("`vname'",16)  _col(54) %8.3f `bin'  " " _col(64) %8.0f   `nc' " " _col(73) %8.0f `nt'
end


********************************************************************************
** Find window increments
********************************************************************************


** NOTE: remove widvar as argument (not used)

capture mata: mata drop findwobs()
mata

void findwobs(real scalar wobs, 
			  real scalar nwin, 
			  real scalar posl, 
			  real scalar posr,
			  runvar,dupsvar,wlengthvar)
{

	st_view(R=.,.,runvar,"`touse'")
	st_view(dups=.,.,dupsvar,"`touse'")
	st_view(Wlength=.,.,wlengthvar,"`touse'")

	N = length(R)
	Nc = sum(R:<0)
	Nt = sum(R:>=0)
	poslold = posl
	posrold = posr

	win = 1

	while((win <= nwin) & (wobs < min((posl,Nt-(posr-Nc-1))))) {

		poslold = posl
		posrold = posr

		while((dups[posl]<wobs) & (sum(R[posl] :<= R[posl..poslold])<wobs)) {
			posl = posl - dups[posl]
		}

		while((dups[posr]<wobs) & (sum(R[posrold..posr] :<= R[posr]))<wobs) {
			posr = posr + dups[posr]
		}

		if (abs(R[posl]) < R[posr]) {
			posl = Nc + 1 - sum(-R[posr] :<= R[1..Nc])
		}
		
		if (abs(R[posl]) > R[posr]) {
			posr = sum(R[(Nc+1)..N] :<= abs(R[posl])) + Nc
		}
		
		wlength = max((-R[posl],R[posr]))
		
		Wlength[posl..poslold] = J(poslold-posl+1,1,wlength)
		Wlength[posrold..posr] = J(posr-posrold+1,1,wlength)

		posl = posl - 1
		posr = posr + 1
		win = win + 1
	}
	
	st_numscalar("posl",posl)
	st_numscalar("posr",posr)
	st_numscalar("wlength",wlength)
}

end


********************************************************************************
** Find recommended length
********************************************************************************

capture mata: mata drop reclength()
mata:
void reclength(Mp,Mw,real scalar level)
{
	if (Mp[1]<level) {
		recl = 0
	} else if (sum(Mp:<level)==0){
		indMw = rows(Mw)
		recl = Mw[indMw]
		ind = rows(Mw)
	} else {
		aux = select(Mw,Mp:<level)
		tmp = min(aux)
		maxindex(Mw:==tmp,1,ind=.,x=.)
		ind = max(ind)
		ind = ind - 1
		recl = Mw[ind]
	}
	st_numscalar("reclength",recl)
	st_numscalar("index",ind)

}
end


********************************************************************************
** Find window length - DEPRECATED: for backward compatibility
********************************************************************************

capture mata: mata drop wlength()
mata:
function wlength(runv, treat, cont, real scalar num)
{
	real scalar xt
	real scalar xc
	
	st_view(R=.,.,runv)
	x0 = min(max(R)\abs(min(R)))
	st_view(X=.,.,runv,treat)
	X = sort(X,1)
	st_view(Y=.,.,runv,cont)
	Z = sort(abs(Y),1)
	m = min(length(X)\length(Z))
	if (m<num) num=m
	xt = X[num]
	xc = Z[num]
	minw = max(xt\xc)
	
	return(xt,xc,minw)
}
end

********************************************************************************
** Find step - DEPRECATED: for backward compatibility
********************************************************************************

capture mata: mata drop findstep()
mata:
void findstep(real scalar minobs,real scalar addobs,real scalar times, runv, treat, cont)
{
	S = wlength(runv,treat,cont,minobs+addobs)-wlength(runv,treat,cont,minobs)
	for(i=1; i<=times-1; i++) {
		U = wlength(runv,treat,cont,minobs+addobs*i)
		L = wlength(runv,treat,cont,minobs+addobs*(i-1))
		Snext = U-L
		S = (S\Snext)
	}
	step = max(S)
	st_numscalar("step",step)
}
end

