
/*******************************************************************************

* Auxiliary program to permute all covariates simultaneously

*!version 0.3 13-Mar-2018

Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare

**** NOTES

* score variable must be recentered at the cutoff
* before running this program:

preserve
keep if `inwindow' `if' `in'

*******************************************************************************/

capture program drop rdwinselect_allcovs
program define rdwinselect_allcovs, rclass
	syntax varlist, treat(string) runvar(string) stat(string) [weights(string)]
	
	if "`weights'"!=""{
		local weight_opt "weights(`weights')"
	}
	
	local nvars: word count `varlist'
	local row = 1
	foreach var of varlist `varlist'{
		rdrandinf_model `var' `treat', stat(`stat') `weight_opt'
		return scalar stat_`row' = r(stat) 
		local ++row
	}
	
end
