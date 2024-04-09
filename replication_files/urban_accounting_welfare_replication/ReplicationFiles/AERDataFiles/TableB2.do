* This program generates TABLE B2 of "Urban Accounting and Welfare" by Klaus Desmet and Esteban Rossi-Hansberg
* October 11, 2012

clear
set memory 12m
set matsize 250

log using TableB2.log, replace
log off


* PARAMETER VALUES FROM MCGRATTAN
global psi=1.4841
global theta=.3358
global totalhours=5110
global delta=0.02
global inter=0.02

* DATASET USED
use DataMSA.dta 

* Total wage income (1000 dollars)
forvalues x=2005/2008 {
	gen totalwages`x'=wage`x'*jobs`x'/1000 
}

* Efficiency wedge
forvalues x=2005/2008 {
	gen efficiencywedge`x'=(gdp`x'/pop`x')/((cap`x'/pop`x')^($theta)*(avhours`x'/$totalhours)^(1-$theta))
}


* Labor productivity wages(1000 dollars per hour worked)
forvalues x=2005/2008 {
	gen laborprodwage`x'=totalwages`x'/agghours`x' 
}

* Labor productivity output(1000 dollars output per hour worked)
forvalues x=2005/2008 {
	gen laborprod`x'=gdp`x'/agghours`x' 
}


* RESHAPE LONG

keep fips efficiencywedge* msaname laborprod* laborprodwage*
reshape long efficiencywedge laborprod laborprodwage, i(fips) j(year)

log on

pwcorr efficiencywedge laborprod, sig
pwcorr efficiencywedge laborprodwage, sig

log close

 

