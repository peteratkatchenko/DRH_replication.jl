* This program generates the inputs needed to run Robustness Exercise 1 for the US, where the elasticity of commuting costs relative to population is 0.25, instead of 0.5
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


* SETS PARAMETER FOR ELASTICITY OF COMMUTING WITH RESPECT TO POPULATION. IN THE STANDARD MODEL GAMMA=1 GIVES US AN ELASTICITY OF 0.5, WITH GAMMA=0.5 IN THIS ROBUSTNESS CHECK WE GET AN ELASTICITY OF 0.5
global gamma=0.5

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

* USE EQUATION (20) TO COMPUTE LOG OF EXCESSIVE FRICTIONS AND ALPHA_5 AND THEN DETERMINE KAPPA
gen loglaborpop=loglabor-(1-$gamma)/2*log(exp(logpop)*1000)
reg loglaborpop dummy2005 dummy2006 dummy2007
scalar alpha5=_b[_cons]
scalar alpha5dum2005=_b[_cons]+_b[dummy2005]
scalar alpha5dum2006=_b[_cons]+_b[dummy2006]
scalar alpha5dum2007=_b[_cons]+_b[dummy2007]

scalar kappa=exp(alpha5-log(2/3)+1/2*log(3.1415))
scalar kappa2005=exp(alpha5dum2005-log(2/3)+1/2*log(3.1415))
scalar kappa2006=exp(alpha5dum2006-log(2/3)+1/2*log(3.1415))
scalar kappa2007=exp(alpha5dum2007-log(2/3)+1/2*log(3.1415))

display kappa
scalar kappa2=round(kappa, 0.00001)

gen logexcfrict=log(exp(loglaborwedge)/kappa2*(2+$gamma)/2*(3.1415/(pop*1000))^($gamma/2))


* KEEP VARIABLES NEEDED TO RUN COUNTERFACTUAL EXERCISE
keep year fips msaname statename logeff logexcfrict pop

reshape wide logeff logexcfrict pop, i(fips) j(year)

drop if pop2005==.
drop if pop2006==.
drop if pop2007==.
drop if pop2008==.







