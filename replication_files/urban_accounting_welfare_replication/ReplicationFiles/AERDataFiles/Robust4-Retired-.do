* This program generates the inputs needed to run Robustness Exercise 4 for the US, where hours worked are computed for people aged 65 and below
* "Urban Accounting and Welfare" by Klaus Desmet and Esteban Rossi-Hansberg
* October 11, 2012

clear
set memory 12m
set matsize 250

* PARAMETER VALUES FROM MCGRATTAN
global psi=1.4841
global theta=.3358
global totalhours=5110
global delta=0.02
global inter=0.02

* DATASET USED
use DataMSA.dta 

* CALCULATE SOME OF THE VARIABLES NEEDED

* only keep 65 and below and change name to standard name

forvalues x=2005/2008 {
	drop avhours`x' avhoursannual70`x'
	rename avhoursannual65`x' avhours`x'
}

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
keep pop* laborwedge* efficiencywedge* medianrent* gdp* fips msaname statename
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

* USE EQUATION (20) TO COMPUTE LOG OF EXCESSIVE FRICTIONS
* Note: we are using kappa estimated in benchmark exercise 0.00176 (InputsNumericalUS.do) 
gen logexcfrict=log(exp(loglaborwedge)/0.00176*3/2*(3.1415/(pop*1000))^0.5)


* KEEP VARIABLES NEEDED TO RUN COUNTERFACTUAL EXERCISE
keep year fips msaname statename logeff logexcfrict pop

reshape wide logeff logexcfrict pop, i(fips) j(year)

drop if pop2005==.
drop if pop2006==.
drop if pop2007==.
drop if pop2008==.
