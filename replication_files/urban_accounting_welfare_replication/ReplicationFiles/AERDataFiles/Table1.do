* This program generates TABLE 1 of "Urban Accounting and Welfare" by Klaus Desmet and Esteban Rossi-Hansberg
* October 11, 2012

clear
set memory 12m
set matsize 250

log using Table1.log, replace

log off

* PARAMETER VALUES FROM MCGRATTAN
global psi=1.4841
global theta=.3358
global totalhours=5110
global delta=0.02
global inter=0.02

* DATASET USED
use DataMSA.dta 

* CALCULATE SOME OF THE VARIABLES NEEDED

* Labor wedge expressed as (1-tau)
forvalues x=2005/2008 {
	gen laborwedge`x'=$psi*privconstot`x'/gdp`x'/(1-avhours`x'/$totalhours)*(avhours`x'/$totalhours)/(1-$theta)
}

* Efficiency wedge
forvalues x=2005/2008 {
	gen efficiencywedge`x'=(gdp`x'/pop`x')/((cap`x'/pop`x')^($theta)*(avhours`x'/$totalhours)^(1-$theta))
}

* DROP IF CERTAIN VARIABLES ARE MISSING
drop if pop2005==.
drop if pop2006==.
drop if pop2007==.
drop if pop2008==.

drop if efficiencywedge2005==.
drop if efficiencywedge2006==.
drop if efficiencywedge2007==.
drop if efficiencywedge2008==.

drop if laborwedge2005==.
drop if laborwedge2006==.
drop if laborwedge2007==.
drop if laborwedge2008==.


* RESHAPE DATA IN LONG FORMAT
keep pop* laborwedge* efficiencywedge* medianrent* gdp* fips
reshape long pop laborwedge efficiencywedge medianrent gdp, i(fips) j(year)

* CREATE TIME DUMMIES
forvalues x=2005/2008 {
	gen dummy`x'=0
	replace dummy`x'=1 if year==`x' 
}

* CREATE LOG VARIABLES 
gen logpop=log(pop)
gen logmedianrent=log(medianrent)
gen logeff=log(efficiencywedge)
* Important: labor wedge is now defined as log of tau 
gen loglaborwedge=log(1-laborwedge)


* DROP VARIABLES THAT ARE MISSING
drop if logpop==.
drop if logeff==.
drop if loglaborwedge==.

log on

* REGRESSION (14) IN PAPER, CORRESPONDING TO BETA_1 IN TABLE 1
reg logpop logeff dummy2005 dummy2006 dummy2007, vce(cluster fips)

log off

* compute predicted population from estimation of regression (14)
scalar coeff1=_b[logeff]
gen predpop=coeff1*logeff
predict respop, residuals

log on
 
* REGRESSION (15) IN PAPER, CORRESPONDING TO BETA_2 IN TABLE 1
* because of bootstrap, standard errors may slightly differ between runs
bootstrap, rep(1000): reg loglaborwedge predpop dummy2005 dummy2006 dummy2007

log off

predict reslaborwedge, residuals
scalar coeff2=_b[predpop]
gen predlaborwedge=coeff2*predpop
egen medpredpop=pctile(predpop), p(50)
gen predlaborwedge2=predlaborwedge-coeff2*medpredpop
gen interreslaborrespop=reslaborwedge*respop

log on

* REGRESSION (16) IN PAPER, INCLUDING INTERACTION TERM, CORRESPONDING TO BETA_3, BETA_4, BETA_5 AND BETA_7 IN TABLE 1
* because of bootstrap, standard errors may slightly differ between runs
reg logmedianrent predlaborwedge2 reslaborwedge respop interreslaborrespop dummy2005 dummy2006 dummy2007, vce(bootstrap, rep(1000))


* REGRESSION (17) IN PAPER, CORRESPONDING TO BETA_6 IN TABLE 1
reg logpop logmedianrent dummy2005 dummy2006 dummy2007, vce(cluster fips) 
test logmedianrent-2=0

log close








