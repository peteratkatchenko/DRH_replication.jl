* This program generates the inputs needed to run the counterfactual exercise for China in "Urban Accounting and Welfare" by Klaus Desmet and Esteban Rossi-Hansberg
* October 11, 2012

clear
set memory 12m
set matsize 250

* PARAMETER VALUES 

global psi=1.5247
global theta=0.5221
global totalhours=5110
global inter=0.2008


* DATASET USED

use DataChina.dta 

rename censpop152005 censpop15
rename censemp2005 censemp
rename emp2005 emp
rename pop2005 pop
rename gdp2005 gdp
rename cons2005 cons


* CALCULATE SOME OF THE VARIABLES NEEDED

* Labor wedge expressed as (1-tau)
gen hours=(censemp*censhours*52)/censpop15/$totalhours
gen laborwedge=$psi/(1-$theta)*cons/gdp*hours/(1-hours)


* Efficiency wedge
gen efficiencywedge=(gdp/pop)^(1-$theta)/($theta/$inter)^$theta/ (hours^(1-$theta))


* CREATE LOG VARIABLES 
gen logpop=log(pop)
gen loglaborwedge=log(1-laborwedge)
gen logeff=log(efficiencywedge)
gen logpop2=log(exp(logpop)*10000)

* USE EQUATION (20) TO COMPUTE LOG OF EXCESSIVE FRICTIONS AND ALPHA_5 AND THEN DETERMINE KAPPA
gen loglaborpop=loglaborwedge-0.5*log(exp(logpop)*10000)
reg loglaborpop 
scalar alpha5=_b[_cons]
scalar kappa=exp(coeffg-ln(2/3)+0.5*ln(3.1415))
display kappa
gen logexcfrict=log(exp(loglaborwedge)/kappa*3/2*(3.1415/(pop*10000))^0.5)

* KEEP VARIABLES NEEDED TO RUN COUNTERFACTUAL EXERCISE
* population in 1,000 to make comparable to U.S. 
replace pop=pop*10
drop if loglaborwedge==.
drop if logpop==.
drop if cityid==.
gen year=2005
keep province city year logeff logexcfrict pop







