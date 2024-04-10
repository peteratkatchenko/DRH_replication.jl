
# This program computes the counter-factual distributions in SA

#Select Country 

using CSV 
using DataFrames
using NLsolve 
using Plots

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
    USBenchmark = CSV.read("C:\\Users\\peter\\.julia\\dev\\dev_econ_replication\\replication_files\\urban_accounting_welfare_replication\\ReplicationFiles\\MatlabPro\\USBenchmark.txt",
    DataFrame) 
    shocks = USBenchmark

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
    ChinaBenchmark = CSV.read("C:\\Users\\peter\\.julia\\dev\\dev_econ_replication\\replication_files\\urban_accounting_welfare_replication\\ReplicationFiles\\MatlabPro\\ChinaBenchmark.txt",
    DataFrame) 
    shocks = ChinaBenchmark

end

Nbar = sum(1000 .*shocks[:,5])
Nlar = max(1000 .*shocks[:,5])

NbarNew = 0

N = 1000 .*shocks[:,5]

c1 = zeros(Float64, 191)
c2 = zeros(Float64, 191)

for i in 1:length(shocks[:,1])

    c1[i] = ubar + (1 + psi)*log((1 + psi))-psi*log(psi)
    c2[i] = (1-theta)*((exp.(shocks[i,2])^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))

    shocks[i,3] = -log(c2[i])+((N[i]*((kappa^2)/pi))^(1/2))*(((1+psi)/c2[i])+(2/3)*exp.(shocks[i,4])) + c1[i]
    shocks[i,3] = log(shocks[i,3])
end 

NbarNew = sum(N)

# Counterfactual efficiency
count = 0

step = 0.1 + China

NbarNewTFP = 0

ubarTFP = ubar

X0 = zeros(Float64, 191)

while abs(Nbar-NbarNewTFP) > 0.02

    for i = 1:length(shocks[:,1])
       
        if China == 0
        
            RTFP = sum(exp.(shocks[:,2])./((shocks[:,5].*1000).^eps).*(shocks[:,5].*1000/Nbar))
            
        else
            
            RTFP = quantile(exp.(shocks[:,2])./((shocks[:,5].*1000).^eps), 0.49)
            
        end    
        
        if NbarNewTFP == 0
            X0 = 1000*shocks[i,5]
        else
            X0 = max(NTFP[i],1)
        end
        
        F=@(R)NsysTFPPA(R)
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off')
        [X,Eval] = fsolve(F,X0,options)
    
        NTFP[i] = max(0, real.(X))
        
        
    end 
 
    NbarNewTFP = sum(NTFP)
    count = count + 1
    if count > 10
        step = step/10
        count = 0
    end
    
    if Nbar-NbarNewTFP > 0.02
        ubarTFP = ubarTFP - step
    else
        ubarTFP = ubarTFP + step
    end
    println(Nbar-NbarNewTFP)
end

NbarNewTFP
max(NTFP)
ubarTFP

# Counterfactual ability
count = 0
step = 0.1 + China
NbarNewA = 0
ubarA = ubar

while abs(Nbar-NbarNewA) > 0.02

    for i = 1:length(shocks[:,1])
        
        if China == 0
        
            AA = sum(exp.(shocks[:,3])./((shocks[:,5].*1000).^ksi).*shocks[:,5].*1000/Nbar)
            
        else
            
            AA = quantile(exp.(shocks[:,3])./((shocks[:,5].*1000).^ksi), 0.49)
            
        end    
        
        if NbarNewA == 0
            X0 = 1000*shocks[i,5]
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

NbarNewA
max(NA)
ubarA

# Counterfactual EF
count = 0
step = 0.1
NbarNewEF = 0
ubarEF = ubar

while abs(Nbar-NbarNewEF) > 0.02

    for i = 1:length(shocks[:,1])
        
        if China == 0
        
            AEF = sum(exp.(shocks[:,4]).*shocks[:,5].*1000/Nbar)
            
        else
            
            AEF = quantile(exp.(shocks[:,4]), 0.49)
            
        end    
        
        if NbarNewEF == 0
            X0 = 1000*shocks[i,5]
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

NbarNewEF
max(NEF)
ubarEF

TMA = sum(max((NA .- N), 0))/Nbar
TMTFP = sum(max((NTFP .- N), 0))/Nbar
TMEF = sum(max((NEF .- N), 0))/Nbar

 
SNR = sort(log.(1000 .*shocks[:,5]))
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


chA = ((NA .- N)./N)
chTFP = ((NTFP .- N)./N)
chEF = ((NEF .- N)./N)
PercChange = hcat(chA, chTFP, chEF)
ShocksModel = shocks[:, 1:5]

if China == 0 
    save("PercChangeUS.jld2", PercChange)
    save("ShocksModelUS.jld2", ShocksModel)
else
    save(PercChangeChina.jld2, PercChange)
    save(ShocksModelChina.jld2, ShocksModel)
end    
    
#Only US
if China == 0
 
    # Counterfactual TFP only
    count = 0
    step = 0.1
    NbarNewTFPO = 0
    ubarTFPO = ubar

    while abs(Nbar-NbarNewTFPO) > 0.02

        for i = 1:length(shocks[:,1])
        
            if NbarNewTFPO == 0
                X0 = 1000*shocks[i,5]
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

    NbarNewTFPO
    max(NTFPO)
    ubarTFPO

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

    while abs(Nbar-NbarNewAO) > 0.02

        for i = 1:length(shocks[:,1])
        
        
            if NbarNewAO == 0
                X0 = 1000*shocks[i,5]
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

    NbarNewAO
    max(NAO)
    ubarAO

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

    while abs(Nbar-NbarNewEFO) > 0.02

        for i = 1:length(shocks[:,1])
        
            if NbarNewEFO == 0
                X0 = 1000*shocks[i,5]
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

    NbarNewEFO
    max(NEFO)
    ubarEFO

    TMEFO = sum(max((NEFO .- N), 0))/Nbar

    SNEFO = sort(log.(NEFO))

    f2p4 = plot(SN, lprob, 
    label="Actual", xlabel="ln(population)", ylabel="ln(prob > population)", legend=:southwest,
    title="'Counterfactual Utility = $(ubarEFO), Reallocation = $(TMEFO)",
    color=:blue, linewidth=2)
    plot!(f2p4, SNEFO, lprob, label="Exc. Frictions Only", legend=:southwest, color=:black, linewidth=2)

    f2 = plot(f2p1, f2p2, f3p3, f4p4, layout=(2,2))


    # Counterfactual Without Shocks
    count = 0
    step = 0.1
    NbarNewNS = 0
    ubarNS = ubar

    while abs(Nbar-NbarNewNS) > 0.02

        for i = 1:length(shocks[:, 1])
        
            if NbarNewNS == 0
                X0 = 1000*shocks[i,5]
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

    NbarNewNS
    max(NNS)
    ubarNS

    SNNS = sort(log.(NNS))

    f3 = plot(SN, lprob, 
    label="Actual", legend=:southwest,
    title="Counterfactual Utility = $(ubarNS)",
    xlabel="ln(population)", ylabel="ln(prob > population)",
    color=:blue, linewidth=2)
    plot(f3, SNNS, lprob, label="No Shocks", legend=:southwest, color=:black, linewidth=2)

end 
