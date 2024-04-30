using StatFiles
using DataFrames 
using Econometrics 
using CSV 
using JLD2
using Statistics


#PARAMETER VALUES FROM MCGRATTAN
psi=1.4841
theta=.3358
totalhours=5110
delta=0.02
inter=0.02

#DATASET USED

filepath2 = joinpath(@__DIR__, "DataMSA.jld2")

data = DataFrame(load(filepath2))

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

data1 = stack(data, Cols(r"pop.*"), [:fips, :msaname, :statename], variable_name=:year, value_name=:pop)
data1.year = parse.(Int, replace.(data1.year, r".*_" => ""))
data1 = sort(data1, [:fips, :year])

data2 = stack(data, Cols(r"laborwedge.*"), [:fips, :msaname, :statename], variable_name=:year, value_name=:laborwedge)
data2.year = parse.(Int, replace.(data2.year, r".*_" => ""))
data2 = sort(data2, [:fips, :year])

data3 = stack(data, Cols(r"efficiencywedge.*"), [:fips, :msaname, :statename], variable_name=:year, value_name=:efficiencywedge)
data3.year = parse.(Int, replace.(data3.year, r".*_" => ""))
data3 = sort(data3, [:fips, :year])

data4 = stack(data, Cols(r"medianrent.*"), [:fips, :msaname, :statename], variable_name=:year, value_name=:medianrent)
data4.year = parse.(Int, replace.(data4.year, r".*_" => ""))
data4 = sort(data4, [:fips, :year])

data5 = stack(data, Cols(r"gdp.*"), [:fips, :msaname, :statename], variable_name=:year, value_name=:gdp)
data5.year = parse.(Int, replace.(data5.year, r".*_" => ""))
data5 = sort(data5, [:fips, :year])

data = hcat(data1, data2, data3, data4, data5, makeunique=true)

data = select(data, [:fips, :msaname, :statename, :year, :pop, :laborwedge, :efficiencywedge, :medianrent, :gdp])

#CREATE TIME DUMMIES

for i in 2005:2008
	data[!, "dummy$i"] = zeros(Int, length(data.year))
	for j in eachindex(data.year)
		if data[j, :year] == i
			data[j, "dummy$i"] = 1
		end 
	end 
end 


#CREATE LOG VARIABLES 

data.logpop = log.(data[!, :pop])

data.logmedianrent = log.(data[!, :medianrent])

data.logeff = log.(data[!, :efficiencywedge])

#Important: labor wedge is now defined as log of tau 

data.loglaborwedge = log.(1 .- data[!, :laborwedge])


#DROP VARIABLES THAT ARE MISSING

dropmissing!(data, [:logpop, :logeff, :loglaborwedge])



#USE EQUATION (20) TO COMPUTE LOG OF EXCESSIVE FRICTIONS AND ALPHA_5 AND THEN DETERMINE KAPPA

data.loglaborpop = data[!, :loglaborwedge] .- 0.5 .* log.(exp.(data[!, :logpop]) .* 1000)

model = fit(EconometricModel, @formula(loglaborpop ~ dummy2005 + dummy2006 + dummy2007), data)

alpha5 = coef(model)[1]

alpha5dum2005 = alpha5 + coef(model)[2]

alpha5dum2006 = alpha5 + coef(model)[3]

alpha5dum2007 = alpha5 + coef(model)[4]

kappa = exp(alpha5 - log(2/3) + (1/2)*log(3.1415))

kappa2005 = exp(alpha5dum2005 - log(2/3) + (1/2)*log(3.1415))

kappa2006 = exp(alpha5dum2006 - log(2/3) + (1/2)*log(3.1415))

kappa2007 = exp(alpha5dum2007 - log(2/3) + (1/2)*log(3.1415))

println(kappa)

kappa2 = round(kappa, digits=5)

data.logexcfrict = log.( exp.(data[!, :loglaborwedge]) ./ kappa2 .* (3/2) .* (3.1415 ./(data[!, :pop] .* 1000)) .^0.5)


#KEEP VARIABLES NEEDED TO RUN COUNTERFACTUAL EXERCISE

data = select(data, [:year, :fips, :msaname, :statename, :logeff, :logexcfrict, :pop])

#RESHAPE DATA FROM LONG TO WIDE

data1 = select(data, [:year, :fips, :msaname, :statename, :logeff])
data1 = unstack(data1, [:fips, :msaname, :statename], :year, :logeff)
rename!(data1, "2005" => "logeff_2005", "2006" => "logeff_2006", "2007" => "logeff_2007", "2008" => "logeff_2008")

data2 = select(data, [:year, :fips, :msaname, :statename, :logexcfrict])
data2 = unstack(data2, [:fips, :msaname, :statename], :year, :logexcfrict)
rename!(data2, "2005" => "logexcfrict_2005", "2006" => "logexcfrict_2006", "2007" => "logexcfrict_2007", "2008" => "logexcfrict_2008")

data3 = select(data, [:year, :fips, :msaname, :statename, :pop])
data3 = unstack(data3, [:fips, :msaname, :statename], :year, :pop)
rename!(data3, "2005" => "pop_2005", "2006" => "pop_2006", "2007" => "pop_2007", "2008" => "pop_2008")

data = hcat(data1, data2, data3, makeunique=true)

data = select(data, :fips, :msaname, :statename, Cols(r"logeff_.*"), Cols(r"logexcfrict_.*"), Cols(r"pop.*"))

dropmissing!(data, Cols(r"pop.*"))

#CALCULATING AVG POP EFF EXCFRICT 

data = select!(data, Cols(r"pop_*"), Cols(r"logeff_*"), Cols(r"logexcfrict_*"))

data.avg_pop = mean.(eachrow(data[:, [:pop_2005, :pop_2006, :pop_2007, :pop_2008]]))

data.avg_logeff = mean.(eachrow(data[:, [:logeff_2005, :logeff_2006, :logeff_2007, :logeff_2008]]))

data.avg_logexcfrict = mean.(eachrow(data[:, [:logexcfrict_2005, :logexcfrict_2006, :logexcfrict_2007, :logexcfrict_2008]]))

data = select!(data, [:avg_pop, :avg_logeff, :avg_logexcfrict])

data.zero = zeros(Float32, 193)

dict = Dict(:avg_pop => data[!, :avg_pop], :avg_logeff => data[!, :avg_logeff], :avg_logexcfrict => data[!, :avg_logexcfrict],
:zero => data[!, :zero])

save("USBenchmark.jld2", dict)