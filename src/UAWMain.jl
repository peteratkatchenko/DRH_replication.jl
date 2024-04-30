module UAWMain 

using CSV 
using DataFrames
using NLsolve 
using Plots

include("NsysTFPPA.jl")
include("NsysAPA.jl")
include("NsysEFPA.jl")
include("NsysTFPOPA.jl")
include("NsysAOPA.jl")
include("NsysEFOPA.jl")
include("NsysNSPA.jl")

import .NsysTFPPA_: NsysTFPPA
import .NsysAPA_: NsysAPA 
import .NsysEFPA_: NsysEFPA 
import .NsysTFPOPA_: NsysTFPOPA 
import .NsysAOPA_: NsysAOPA 
import .NsysEFOPA_: NsysEFOPA 
import .NsysNSPA_: NsysNSPA 

#Select Country 
China = 0 #Make 0 if US and 1 if China

if China == 0
    
    #US Parameters
    psi = 1.4841
    theta = .3358
    totalhours = 5110
    beta = 0.98
    delta = 0.02
    inter = 0.02
    eps = .00 #(called omega in main text)
    ksi = .00
    ubar = 10
    kappa = 0.002
    USBenchmark = CSV.read("C:\\Users\\peter\\.julia\\dev\\USBenchmark.csv",
    DataFrame) 
    shocks = USBenchmark

    dictmain = Dict(:psi => psi, :theta => theta, :totalhours => totalhours, :beta => beta,
    :delta => delta, :inter => inter, :eps => eps, :ksi => ksi, :ubar => ubar, :kappa => kappa,
    :shocks => shocks)
else
    
    #China Parameters
    psi = 1.52
    theta = .52
    totalhours = 5110
    beta = 0.90
    delta = 0.2
    inter = 0.2
    eps = .00
    ksi = .00
    ubar = 10
    kappa = 0.001
    ChinaBenchmark = CSV.read("C:\\Users\\peter\\.julia\\dev\\dev_econ_replication\\replication_files\\urban_accounting_welfare_replication\\ReplicationFiles\\MatlabPro\\ChinaBenchmark.csv",
    DataFrame) 
    shocks = ChinaBenchmark

    dictmain = Dict(:psi => psi, :theta => theta, :totalhours => totalhours, :beta => beta,
    :delta => delta, :inter => inter, :eps => eps, :ksi => ksi, :ubar => ubar, :kappa => kappa,
    :shocks => shocks)
end

psi = dictmain[:psi]
theta = dictmain[:theta]
totalhours = dictmain[:totalhours]
beta = dictmain[:beta]
delta = dictmain[:delta]
inter = dictmain[:inter]
eps = dictmain[:eps]
ksi = dictmain[:ksi]
ubar = dictmain[:ubar]
kappa = dictmain[:kappa]
shocks = dictmain[:shocks]

Nbar = sum(1000 .*shocks[:,1])
Nlar = maximum(1000 .*shocks[:,1])

N = 1000 .*shocks[:,1]

c1 = zeros(Float64, 193)
c2 = zeros(Float64, 193)

for i in 1:length(shocks[:,1])

    c1[i] = ubar + (1 + psi)*log((1 + psi))-psi*log(psi)
    c2[i] = (1-theta)*((exp(shocks[i,2])^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))

    shocks[i,4] = -log(c2[i])+((N[i]*((kappa^2)/pi))^(1/2))*(((1+psi)/c2[i])+(2/3)*exp(shocks[i,3])) + c1[i]
    shocks[i,4] = log(shocks[i,4])
end 

NbarNew::Float64 = sum(N)

# Counterfactual efficiency
count = 0

step = 0.1 + China

NbarNewTFP = 0

ubarTFP = ubar

NTFP = zeros(Float64, 193)

dictCE = Dict()

dictCE[:ubarTFP] = ubarTFP


while abs(Nbar-NbarNewTFP) > 0.02

    for i = 1:length(shocks[:,1])

        D = zeros(Float64, 1)

        X0 = zeros(Float64, 1)

        if China == 0
        
            RTFP = sum(exp.(shocks[:,2])./((shocks[:,1].*1000).^eps).*(shocks[:,1].*1000 ./Nbar))
            dictCE[:RTFP] = RTFP
        else
            
            RTFP = quantile(exp.(shocks[:,2])./((shocks[:,1].*1000).^eps), 0.49)
            dictCE[:RTFP] = RTFP
        end    

        
        if NbarNewTFP == 0
            X0[1] = 1000*shocks[i,1]
        else
            X0[1] = max(NTFP[i], 1)
        end

        objfun(D, X) = NsysTFPPA(D, X, i, dictmain, dictCE)
        
        sol = nlsolve(objfun, X0, iterations=1000000, ftol=0.0001)
        X = sol.zero
        X_ = X[1]
    
        NTFP[i] = max(0, X_)
        
        
    end 
 
    global NbarNewTFP = sum(NTFP)
    global count = count + 1
    if count > 10
       global step = step/10
       global count = 0
    end
    
    if Nbar-NbarNewTFP > 0.02
        global ubarTFP = ubarTFP - step
    else
        global ubarTFP = ubarTFP + step
    end
    println(Nbar-NbarNewTFP)
end

println(NbarNewTFP, max(NTFP), ubarTFP)

#=

# Counterfactual ability
count = 0

step = 0.1 + China

NbarNewA = 0

ubarA = ubar

NA = zeros(Float64, 193)

while abs(Nbar-NbarNewA) > 0.02

    for i = 1:length(shocks[:,1])
        
        if China == 0
        
            AA = sum(exp.(shocks[:,4])./((shocks[:,1].*1000).^ksi).*shocks[:,1].*1000 ./Nbar)
            
        else
            
            AA = quantile(exp.(shocks[:,4])./((shocks[:,1].*1000).^ksi), 0.49)
            
        end    
        
        if NbarNewA == 0
            X0 = 1000*shocks[i,1]
        else
            X0 = max(NA[i], 1)
        end
        
        F=@(R)NsysAPA(R)
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off')
        [X,Eval] = fsolve(F, X0, options)
    
        NA[i] = max(0, real.(X))
            
    end 
 
    NbarNewA = sum(NA)
    count = count + 1
    if count > 10
        step = step/10
        count = 0
    end
    
    if Nbar-NbarNewA > 0.02
        ubarA = ubarA - step
    else
        ubarA = ubarA + step
    end
    println(Nbar-NbarNewA)
end

println(NbarNewA, max(NA), ubarA)

# Counterfactual EF
count = 0

step = 0.1

NbarNewEF = 0

ubarEF = ubar

NEF = zeros(Float64, 193)

while abs(Nbar-NbarNewEF) > 0.02

    for i = 1:length(shocks[:,1])
        
        if China == 0
        
            AEF = sum(exp.(shocks[:,3]).*shocks[:,1].*1000 ./Nbar)
            
        else
            
            AEF = quantile(exp.(shocks[:,3]), 0.49)
            
        end    
        
        if NbarNewEF == 0
            X0 = 1000*shocks[i,1]
        else
            X0 = max(NEF[i], 1)
        end
        
        F=@(R)NsysEFPA(R)
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off')
        [X,Eval] = fsolve(F,X0,options)
    
        NEF[i] = max(0, real.(X))
        
    end 
 
    NbarNewEF = sum(NEF)
    count = count + 1
    if count > 10
        step = step/10
        count = 0
    end
    
    if Nbar-NbarNewEF > 0.02
        ubarEF = ubarEF - step
    else
        ubarEF = ubarEF + step
    end
     println(Nbar-NbarNewEF)
end

println(NbarNewEF, max(NEF), ubarEF)

TMA = sum(max((NA .- N), 0))/Nbar
TMTFP = sum(max((NTFP .- N), 0))/Nbar
TMEF = sum(max((NEF .- N), 0))/Nbar

 
SNR = sort(log.(1000 .*shocks[:,1]))
SN = sort(log.(N))
SNTFP = sort(log.(NTFP))
SNA = sort(log.(NA))
SNEF = sort(log.(NEF))

lprob = zeros(Float64, length(N))

for i = 1:length(N)
    lprob[i] = log((length(N) - i + 1)/length(N))
end

f1p1 = plot(SNR, lprob, color=:red, linewidth=2, label="Actual", legend=:southwest, title="Model Utility = $(ubar)",
xlabel="ln(population)", ylabel="ln(prob > population)")
plot!(p1, SN, lprob, color=:blue, linewidth=2, label="Modeled", legend=:southwest)

f1p2 = plot(SN, lprob, label="Actual", legend=:southwest, 
title="Counterfactual Utility = $(ubarTFP), Reallocation = $(TMTFP)", 
xlabel="ln(population)",
ylabel="ln(prob > population)",
color=:blue, linewidth=2)
plot!(p2, SNTFP, lprob, label="Avg. Efficiency", legend=:southwest, color=:green, linewidth=2)

f1p3 = plot(SN, lprob, label="Actual", legend=:southwest,
title="Counterfactual Utility = $(ubarA), Reallocation = $(TMA)",
xlabel="ln(population)", ylabel="ln(prob > population)",
color=:blue, linewidth=2)
plot!(p3, SNA, lprob, label="Avg. Amenities", legend=:southwest, color=:magenta, linewidth=2)

f1p4 = plot(SN,lprob, 
label="Actual", legend=:southwest,
xlabel="ln(population)",
ylabel="ln(prob > population)",
title="Counterfactual Utility = ',num2str(ubarEF), Reallocation = $(TMEF)",
annotation=(15.5, 1, "Counterfactuals Without One Shock, κ = $(kappa), ω = $(eps), ζ = $(ksi)", :left),
color=:blue, linewidth=2)
plot(p4, SNEF, lprob, label="Exc. Frictions", legend=:southwest, color=:black, linewidth=2)


f1 = plot(f1p1, f1p2, f1p3, f1p4, layout=(2,2))
display(f1)

chA = ((NA .- N)./N)
chTFP = ((NTFP .- N)./N)
chEF = ((NEF .- N)./N)
PercChange = hcat(chA, chTFP, chEF)
ShocksModel = shocks[:, 1:4]

if China == 0 
    CSV.write("PercChangeUS.csv", PercChange)
    CSV.write("ShocksModelUS.csv", ShocksModel)
else
    save("PercChangeChina.csv", PercChange)
    save("ShocksModelChina.csv", ShocksModel)
end    
    
#Only US
if China == 0
 
    # Counterfactual TFP only
    count = 0
    step = 0.1
    NbarNewTFPO = 0
    ubarTFPO = ubar

    NTFPO = zeros(Float64, 193)

    while abs(Nbar-NbarNewTFPO) > 0.02

        for i = 1:length(shocks[:,1])
        
            if NbarNewTFPO == 0
                X0 = 1000*shocks[i,1]
            else
                X0 = max(NTFPO[i], 1)
            end
        
            F=@(R)NsysTFPOPA(R)
            options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off')
            [X,Eval] = fsolve(F,X0,options)
    
            NTFPO[i] = max(0, real.(X))
        
        end 
 
        NbarNewTFPO = sum(NTFPO)
        count = count + 1
        if count > 10
            step = step/10
            count = 0
        end
    
        if Nbar-NbarNewTFPO > 0.02
            ubarTFPO = ubarTFPO - step
        else
            ubarTFPO = ubarTFPO + step
        end
        println(Nbar-NbarNewTFPO)
    end

    println(NbarNewTFPO, max(NTFPO), ubarTFPO)

    TMTFPO = sum(max((NTFPO .- N), 0))/Nbar

   
    f2p1 = plot(SNR, lprob, label="Actual", legend=:southwest,
    xlabel="ln(population)", ylabel="ln(prob > population)",
    title="Model Utility = $(ubar)", 
    annotation=(15.5, 1, "Counterfactuals with Only One Shock, κ = $(kappa), ω = $(eps), ζ = $(ksi)"),
    color=:red, linewidth=2)
    plot!(f2p1, SN, lprob, label="Modeled", legend=:southwest, color=:blue, linewidth=2)

    SNTFPO = sort(log.(NTFPO))

    f2p2 = plot(SN, lprob,
    xlabel="ln(population)", ylabel="ln(prob > population)",
    label="Actual",
    legend=:southwest,
    title="Counterfactual Utility = $(ubarTFPO), Reallocation = $(TMTFPO)",
    color=:blue, linewidth=2)
    plot(f2p2, SNTFPO, lprob, label="Efficiency Only", legend=:southwest, color=:green, linewidth=2)

    # Counterfactual Amenities Only
    count = 0
    step = 0.1
    NbarNewAO = 0
    ubarAO = ubar

    NAO = zeros(Float64, 193)

    while abs(Nbar-NbarNewAO) > 0.02

        for i = 1:length(shocks[:,1])
        
        
            if NbarNewAO == 0
                X0 = 1000*shocks[i,1]
            else
                X0 = max(NAO[i], 1)
            end
        
            F=@(R)NsysAOPA(R)
            options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off')
            [X,Eval] = fsolve(F,X0,options)
    
            NAO[i] = max(0, real.(X))
        
        end 
 
        NbarNewAO = sum(NAO)
        count = count + 1
        if count > 10
            step = step/10
            count = 0
        end
    
        if Nbar-NbarNewAO > 0.02
            ubarAO = ubarAO - step
        else
            ubarAO = ubarAO + step
        end
        println(Nbar-NbarNewAO)
    end

    println(NbarNewAO, max(NAO), ubarAO)

    TMAO = sum(max((NAO .- N), 0))/Nbar

    SNAO = sort(log.(NAO))

    f2p3 = plot(SN, lprob, 
    xlabel="ln(population)", ylabel="ln(prob > population)",
    legend=:southwest,
    label="Actual",
    title="Counterfactual Utility = $(ubarAO), Reallocation = $(TMAO)",
    color=:blue, linewidth=2)
    plot!(f2p3, SNAO, lprob, label="Amenities Only", legend=:southwest, color=:magenta, linewidth=2)

    # Counterfactual Excessive Frictions Only
    count = 0
    step = 0.1
    NbarNewEFO = 0
    ubarEFO = ubar

    NEFO = zeros(Float64, 193)

    while abs(Nbar-NbarNewEFO) > 0.02

        for i = 1:length(shocks[:,1])
        
            if NbarNewEFO == 0
                X0 = 1000*shocks[i,1]
            else
                X0 = max(NEFO[i], 1)
            end
        
            F=@(R)NsysEFOPA(R)
            options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off')
            [X,Eval] = fsolve(F,X0,options)
    
            NEFO[i] = max(0, real.(X))
        
        end 
 
        NbarNewEFO = sum(NEFO)
        count = count + 1
        if count > 10
            step = step/10
            count = 0
        end
    
        if Nbar-NbarNewEFO > 0.02
            ubarEFO = ubarEFO - step
        else
            ubarEFO = ubarEFO + step
        end
        println(Nbar-NbarNewEFO)
    end

    println(NbarNewEFO, max(NEFO), ubarEFO)

    TMEFO = sum(max((NEFO .- N), 0))/Nbar

    SNEFO = sort(log.(NEFO))

    f2p4 = plot(SN, lprob, 
    label="Actual", xlabel="ln(population)", ylabel="ln(prob > population)", legend=:southwest,
    title="'Counterfactual Utility = $(ubarEFO), Reallocation = $(TMEFO)",
    color=:blue, linewidth=2)
    plot!(f2p4, SNEFO, lprob, label="Exc. Frictions Only", legend=:southwest, color=:black, linewidth=2)

    f2 = plot(f2p1, f2p2, f3p3, f4p4, layout=(2,2))
    display(f2)

    # Counterfactual Without Shocks
    count = 0
    step = 0.1
    NbarNewNS = 0
    ubarNS = ubar

    NNS = zeros(Float64, 193)

    while abs(Nbar-NbarNewNS) > 0.02

        for i = 1:length(shocks[:, 1])
        
            if NbarNewNS == 0
                X0 = 1000*shocks[i,1]
            else
                X0 = max(NNS[i], 1)
            end
        
            F=@(R)NsysNSPA(R)
            options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off')
            [X,Eval] = fsolve(F,X0,options)
    
            NNS[i] = max(0, real.(X))
        
        end 
 
        NbarNewNS = sum(NNS)
        count = count + 1
        if count > 10
            step = step/10
            count = 0
        end
    
        if Nbar-NbarNewNS > 0.02
            ubarNS = ubarNS - step
        else
            ubarNS = ubarNS + step
        end
        println(Nbar-NbarNewNS)
    end

    println(NbarNewNS, max(NNS), ubarNS)

    SNNS = sort(log.(NNS))

    f3 = plot(SN, lprob, 
    label="Actual", legend=:southwest,
    title="Counterfactual Utility = $(ubarNS)",
    xlabel="ln(population)", ylabel="ln(prob > population)",
    color=:blue, linewidth=2)
    plot(f3, SNNS, lprob, label="No Shocks", legend=:southwest, color=:black, linewidth=2)
    display(f3)
end 
=#


end #End of module 