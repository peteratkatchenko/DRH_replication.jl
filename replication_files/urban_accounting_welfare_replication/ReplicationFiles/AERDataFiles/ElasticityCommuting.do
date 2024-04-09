* This program generates the elasticity of commuting costs relative to population, based on incorporated places with a populaton of more than 100,000
* "Urban Accounting and Welfare" by Klaus Desmet and Esteban Rossi-Hansberg
* October 11, 2012

clear
set more off
set memory 512m

log using elasticitycommuting.log, replace

log off

* INCORPORATED PLACES WITH POPULATION >= 100,000
use elasticitycities.dta

drop if pop < 100000
sum pop
more

gen logpop=log(population)
gen logarea=log(arealand)

reg logarea logpop
scalar coeff1=round(_b[logpop],0.1)
scalar elasticitycommplaces=coeff1/2

log on

* Elasticity of population to commuting, based on incorporated places with population >= 100,000
display elasticitycommplaces

log off

clear

* MSA USED IN NUMERICAL PART OF PAPER
use elasticitymsa.dta

gen logpop=log(pop)
gen logarea=log(area)

reg logarea logpop
scalar coeff1=round(_b[logpop],0.1)
scalar elasticitycommmsa=coeff1/2

log on

* Elasticity of population to commuting, based on MSAs used in paper
display elasticitycommmsa

log close
