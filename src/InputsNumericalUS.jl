using StatFiles
using DataFrames 
using Econometrics 
using CSV 
using RCall

#PARAMETER VALUES FROM MCGRATTAN
psi=1.4841
theta=.3358
totalhours=5110
delta=0.02
inter=0.02

#DATASET USED
data =  DataFrame(load("C:\\Users\\peter\\.julia\\dev\\dev_econ_replication\\replication_files\\urban_accounting_welfare_replication\\ReplicationFiles\\AERDataFiles\\DataMSA.dta")) 

#Changing variable names to facilitate column splitting
rename!(data, "pop2005" => "pop_2005", "pop2006" => "pop_2006", "pop2007" => "pop_2007", "pop2008" => "pop_2008")
rename!(data, "medianrent2005" => "medianrent_2005", "medianrent2006" => "medianrent_2006", "medianrent2007" => "medianrent_2007", "medianrent2008" => "medianrent_2008")
rename!(data, "gdp2005" => "gdp_2005", "gdp2006" => "gdp_2006", "gdp2007" => "gdp_2007", "gdp2008" => "gdp_2008")


#data[175, :cap2005] and data[220, :cap2005] should be positive
data[175, :cap2005] = -(data[175, :cap2005])

data[220, :cap2005] = -(data[220, :cap2005])

#Missing string values in :msaname are stored as empty spaces 
replace!(data[!, :msaname], "" => missing)

#CALCULATE SOME OF THE VARIABLES NEEDED

#Labor wedge expressed as (1-tau)
for i in 2005:2008  
	data[!, Symbol("laborwedge_$i")] = psi .* data[!, Symbol("privconstot$i")] ./ data[!, Symbol("gdp_$i")] ./ (1 .- data[!, Symbol("avhours$i")] ./ totalhours) .* (data[!, Symbol("avhours$i")] ./ totalhours) ./ (1-theta)
end

#Efficiency wedge
for i in 2005:2008 
	data[!, Symbol("efficiencywedge_$i")] = (data[!, Symbol("gdp_$i")] ./ data[!, Symbol("pop_$i")]) ./ (((data[!, Symbol("cap$i")] ./ data[!, Symbol("pop_$i")]) .^(theta)) .* ((data[!, Symbol("avhours$i")] ./totalhours) .^(1-theta)))
end 

#DROP IF CERTAIN VARIABLES ARE MISSING

dropmissing!(data, Cols(r"pop.*"))

dropmissing!(data, Cols(r"efficiencywedge.*"))

dropmissing!(data, Cols(r"laborwedge.*"))

#RESHAPE DATA IN LONG FORMAT

data = select(data, r"pop.*", r"laborwedge.*", r"efficiencywedge.*", r"medianrent.*", r"gdp.*", :fips, :msaname, :statename)




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