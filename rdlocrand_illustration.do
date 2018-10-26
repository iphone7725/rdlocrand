********************************************************************************
** RDLOCRAND Stata Package 
** Do-file for Empirical Illustration
** Authors: Matias D. Cattaneo, Rocio Titiunik and Gonzalo Vazquez-Bare
** Last update: 13-Mar-2018
********************************************************************************
** hlp2winpdf, cdn(rdrandinf) replace
** hlp2winpdf, cdn(rdwinselect) replace
** hlp2winpdf, cdn(rdsensitivity) replace
** hlp2winpdf, cdn(rdrbounds) replace
********************************************************************************
clear all
set more off
set linesize 80

**************************************************************************
** Summary Stats
**************************************************************************
use rdlocrand_senate.dta, clear
global covariates presdemvoteshlag1 population demvoteshlag1 ///
                  demvoteshlag2 demwinprv1 demwinprv2 dopen dmidterm
describe $covariates
summarize demmv $covariates

gen D = demmv>=0
ttest demvoteshfor2, by(D)

timer clear

**************************************************************************
** rdwinselect
**************************************************************************
** Window selection with default options
timer on 1
rdwinselect demmv $covariates, cutoff(0)
timer off 1
timer list 1

** Window selection setting window length and increments (replicate CFT)
timer on 2
rdwinselect demmv $covariates, wmin(.5) wstep(.125) reps(10000)
timer off 2
timer list 2

** Window selection using large sample approximation and plotting p-values
quietly rdwinselect demmv $covariates, wmin(.5) wstep(.125) ///
                    nwin(80) approximate plot

**************************************************************************
** rdrandinf
**************************************************************************
** Randomization inference using recommended window
timer on 4
rdrandinf demvoteshfor2 demmv, wl(-.75) wr(.75)
timer off 4
timer list 4

** Randomization inference using recommended window, all statistics
timer on 5
rdrandinf demvoteshfor2 demmv, wl(-.75) wr(.75) statistic(all)
timer off 5
timer list 5

** Randomization inference using recommended window using rdwinselect
timer on 6
rdrandinf demvoteshfor2 demmv, statistic(all) covariates($covariates) /*
	*/ wmin(.5) wstep(.125) level(0.16) quietly rdwreps(10000)
timer off 6
timer list 6

** Randomization inference using recommended window, linear adjustment
timer on 7
rdrandinf demvoteshfor2 demmv, statistic(all) wl(-.75) wr(.75) p(1)
timer off 7
timer list 7

** Randomization inference under interference
timer on 8
rdrandinf demvoteshfor2 demmv, wl(-.75) wr(.75) interfci(.05)
timer off 8
timer list 8

**************************************************************************
** rdsensitivity
**************************************************************************
timer on 9
rdsensitivity demvoteshfor2 demmv, wlist(.75(.25)2) tlist(0(1)20) nodots verbose
timer off 9
timer list 9

** Obtain 95 percent confidence interval for window [-.75 ; .75]
timer on 10
rdsensitivity demvoteshfor2 demmv, wlist(.75(.25)2) tlist(0(1)20) nodots ci(.75)
timer off 10
timer list 10

** Replicate contour plot
timer on 11
rdsensitivity demvoteshfor2 demmv, wlist(.75(.25)10) tlist(0(1)20) nodots ///
                                   saving(output/graphdata)
preserve
use output/graphdata, clear
twoway contour pvalue t w, ccuts(0(0.05)1)
restore
timer off 11
timer list 11

preserve
use output/graphdata, clear
twoway contour pvalue t w, ccuts(0(0.05)1) ccolors(gray*0.01 gray*0.05 /*
	*/ gray*0.1 gray*0.15 gray*0.2 gray*0.25 gray*0.3 gray*0.35 /*
	*/ gray*0.4 gray*0.5 gray*0.6 gray*0.7 gray*0.8 gray*0.9 gray /*
	*/ black*0.5  black*0.6 black*0.7 black*0.8 black*0.9 black) /*
	*/ xlabel(.75(1.25)10) ylabel(0(2)20, nogrid) graphregion(fcolor(none))
restore

** rdsensitivity to calculate CI from within rdrandinf
timer on 13
rdrandinf demvoteshfor2 demmv, wl(-.75) wr(.75) ci(.05 3(1)20)
timer off 13
timer list 13

** VERY SLOW:
** rdsensitivity demvoteshfor2 demmv, wlist(.75) tlist(3(1)20) ///
**	                                  nodraw nodots ci(.75)


**************************************************************************
** rdrbounds
**************************************************************************
timer on 14
rdrbounds demvoteshfor2 demmv, expgamma(1.5 2 3) wlist(.5 .75 1) reps(1000)
timer off 14
timer list 14

** Bernoulli and fixed margins p-values
timer on 15
rdrbounds demvoteshfor2 demmv, expgamma(1.5 2 3) wlist(.5 .75 1) reps(1000) fmpval
timer off 15
timer list 15

**************************************************************************
** rdrandinf with eval options
**************************************************************************
qui sum demmv if abs(demmv)<=.75 & demmv>=0 & demmv!=. & demvoteshfor2!=.
local mt=r(mean)
qui sum demmv if abs(demmv)<=.75 & demmv<0  & demmv!=. & demvoteshfor2!=.
local mc=r(mean)
rdrandinf demvoteshfor2 demmv, wl(-.75) wr(.75) p(1) evall(`mc') evalr(`mt')


