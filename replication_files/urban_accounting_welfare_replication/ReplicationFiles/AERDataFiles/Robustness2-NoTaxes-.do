* This program generates the inputs needed to run Robustness Exercise 2 for the US, where the labor wedge ignores the part due to taxes
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

* Consumption tax 
forvalues x=2005/2006 {
	gen constaxtot`x'=salestax`x'/privconstot`x' 
}
gen constaxtot2007=constaxtot2006
gen constaxtot2008=constaxtot2006

* Labor wedge expressed as (1-tau). It is computed NET OF TAXES, meaning that we compute THAT PART OF THE LABOR WEDGE THAT CANNOT BE EXPLAINED BY VARIATIONS IN LOCAL TAXES
forvalues x=2005/2008 {
	gen laborwedge`x'=$psi*(1+constaxtot`x')*privconstot`x'/gdp`x'/(1-avhours`x'/$totalhours)*(avhours`x'/$totalhours)/(1-$theta)*(1/(1-inctax`x'))
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
gen loglaborpop=loglaborwedge-0.5*log(exp(logpop)*1000)
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

gen logexcfrict=log(exp(loglaborwedge)/kappa2*3/2*(3.1415/(pop*1000))^0.5)


* KEEP VARIABLES NEEDED TO RUN COUNTERFACTUAL EXERCISE
keep year fips msaname statename logeff logexcfrict pop

reshape wide logeff logexcfrict pop, i(fips) j(year)

drop if pop2005==.
drop if pop2006==.
drop if pop2007==.
drop if pop2008==.







