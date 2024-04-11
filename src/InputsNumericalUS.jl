using StatFiles
using DataFrames 
using Econometrics 
using CSV 

#PARAMETER VALUES FROM MCGRATTAN
psi=1.4841
theta=.3358
totalhours=5110
delta=0.02
inter=0.02

#DATASET USED
data =  DataFrame(load("C:\\Users\\peter\\.julia\\dev\\dev_econ_replication\\replication_files\\urban_accounting_welfare_replication\\ReplicationFiles\\AERDataFiles\\DataMSA.dta")) 

#data[175, :cap2005] and data[220, :cap2005] should be positive - changing to positive numbers

#CALCULATE SOME OF THE VARIABLES NEEDED

#Labor wedge expressed as (1-tau)
for i in 2005:2008  
	data[!, Symbol("laborwedge$i")] = psi .* data[!, Symbol("privconstot$i")] ./ data[!, Symbol("gdp$i")] ./ (1 .- data[!, Symbol("avhours$i")] ./ totalhours) .* (data[!, Symbol("avhours$i")] ./ totalhours) ./ (1-theta)
end

#Efficiency wedge
for i in 2005:2008 
	data[!, Symbol("efficiencywedge$i")] = #=(data[!, Symbol("gdp$i")] ./ data[!, Symbol("pop$i")]) ./ ((=#(data[!, Symbol("cap$i")] ./ data[!, Symbol("pop$i")]) #=.^(theta)) .* ((data[!, Symbol("avhours$i")] ./totalhours) .^(1-theta)))=#
end 


CSV.write("data.csv", data)

for i in 1:nrow(data)
    if !ismissing(data[i, :efficiencywedge2005]) == true  
        if data[i, :efficiencywedge2005] < 0
            return println("Negative value")
        end 
        println("No negative values")
    end    
end 



#DROP IF CERTAIN VARIABLES ARE MISSING
filter!(data, row -> ismissing(row.pop2005) || 
ismissing(row.pop2006) ||
ismissing(row.pop2007) ||
ismissing(row.pop2008))

filter!(data, row -> ismissing(row.efficiencywedge2005) ||
ismissing(efficiencywedge2006) ||
ismissing(efficiencywedge2007) ||
ismissing(efficiencywedge2008))

filter!(data, row -> ismissing(row.laborwedge2005) ||
ismissing(laborwedge2006) ||
ismissing(laborwedge2007) ||
ismissing(laborwedge2008))


#RESHAPE DATA IN LONG FORMAT

data = select(data, r"pop.*" || r"laborwedge.*" || r"medianrent.*" || r"gdp.*", :fips, :msaname, :statename)


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


log close








