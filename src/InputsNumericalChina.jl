using StatFiles
using DataFrames 
using Econometrics 
using CSV 

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
"censhours2005" => "censhours",
"emp2005" => "emp",
"pop2005" => "pop",
"gdp2005" => "gdp",
"cons2005" => "cons")

#CALCULATE SOME OF THE VARIABLES NEEDED

#Labor wedge expressed as (1-tau)

function hours_fun(censemp::Vector, censhours::Vector, censpop15::Vector)
(censemp .*censhours .*52) ./(censpop15 ./ totalhours)
end 
 
data = transform(data, [:censemp, :censhours, :censpop15] => hours_fun => :hours) 

function laborwedge_fun(cons::Vector, gdp::Vector, hours::Vector)
(psi/(1-theta) .* cons) ./(gdp .* hours)./(1 .- hours)
end

data = transform(data, [:cons, :gdp, :hours] => laborwedge_fun => :laborwedge)
    

#Efficiency wedge

function effwedge_fun(gdp::Vector, pop::Vector, hours::Vector)
(gdp ./pop) .^(1-theta) ./(theta/inter)^theta ./ (hours .^(1-theta))
end 

data = transform(data,  [:gdp, :pop, :hours] => effwedge_fun => :efficiencywedge)
  
 

#CREATE LOG VARIABLES

data.logpop = log.(data[!, :pop]) 

data.loglaborwedge = log.(1 .- data[!, :laborwedge])
 
data.logeff = log.(data[!, :efficiencywedge]) 

data.logpop2 = log.(exp.(data[!, :logpop]).*10000) 


#USE EQUATION (20) TO COMPUTE LOG OF EXCESSIVE FRICTIONS AND ALPHA_5 AND THEN DETERMINE KAPPA

data.loglaborpop = data[!, :loglaborwedge] .- 0.5 .* log.(exp.(data[!, :logpop]) .*10000)
 
model = fit(EconometricModel, @formula(loglaborpop ~ PLACEHOLDER), data)

alpha5 = coef(model)[1] 

kappa = exp(coeffg - ln(2/3)+0.5*ln(3.1415))

println(kappa)

data.logexcfrict = log.(exp.(data[!, :loglaborwedge]) ./ kappa .* (3/2) .* (3.1415 ./ (data[!, :pop] .* 10000)).^0.5)

#KEEP VARIABLES NEEDED TO RUN COUNTERFACTUAL EXERCISE
#population in 1,000 to make comparable to U.S. 

data.pop = data[!, :pop] .* 10

filter!(row -> ismissing(row.loglaborwedge) ,data)

filter!(row -> ismissing(row.logpop), data)

filter!(row -> ismissing(row.cityid), data)

data.year = fill(2005, nrow(data))

data = select(data, [:province, :city, :year, :logeff, :logexcfrict, :pop])