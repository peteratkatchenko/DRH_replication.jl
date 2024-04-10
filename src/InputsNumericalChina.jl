using StatFiles
using DataFrames 

#PARAMETER VALUES 

psi=1.5247
theta=0.5221
totalhours=5110
inter=0.2008

#DATASET USED

data = DataFrame(load("C:\\Users\\peter\\.julia\\dev\\dev_econ_replication\\replication_files\\urban_accounting_welfare_replication\\ReplicationFiles\\AERDataFiles\\DataChina.dta")) 

rename!(data, 
"censpop152005" => "censpop15",
"censemp2005" => "censemp",
"emp2005" => "emp",
"pop2005" => "pop",
"gdp2005" => "gdp",
"cons2005" => "cons")


#CALCULATE SOME OF THE VARIABLES NEEDED

#Labor wedge expressed as (1-tau)
transform(data, (censemp .*censhours .*52) ./(censpop15 ./totalhours))

gen hours
gen laborwedge=$psi/(1-$theta)*cons/gdp*hours/(1-hours)


#Efficiency wedge
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







