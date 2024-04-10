using StatFiles
using DataFrames 
using Econometrics 

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

function hours_fun(x::Vector, y::Vector, z::Vector)
(x .*y .*52) ./(z ./ totalhours)
end 

transform(data, [:censemp, :censhours, :censpop15] => hours_fun => :hours)


transform(data, psi/(1-theta)*cons/gdp*hours/(1-hours) => laborwedge)



#Efficiency wedge
transform(data, (gdp/pop)^(1-$theta)/($theta/$inter)^$theta/ (hours^(1-$theta)) => efficiencywedge)



#CREATE LOG VARIABLES
transform(data, :pop => log => :logpop)
transform(data, (1 .- laborwedge) => log => loglaborwedge)
transform(data, efficiencywedge => log => logeff)
transform(data, exp(:logpop)*10000 => log => logpop2)


#USE EQUATION (20) TO COMPUTE LOG OF EXCESSIVE FRICTIONS AND ALPHA_5 AND THEN DETERMINE KAPPA

gen loglaborpop=loglaborwedge-0.5*log(exp(logpop)*10000)

reg loglaborpop 
scalar alpha5=_b[_cons]
scalar kappa=exp(coeffg-ln(2/3)+0.5*ln(3.1415))
display kappa
gen logexcfrict=log(exp(loglaborwedge)/kappa*3/2*(3.1415/(pop*10000))^0.5)

#KEEP VARIABLES NEEDED TO RUN COUNTERFACTUAL EXERCISE
#population in 1,000 to make comparable to U.S. 

replace pop=pop*10
drop if loglaborwedge==.
drop if logpop==.
drop if cityid==.
gen year=2005
keep province city year logeff logexcfrict pop







