using CSV 
using DataFrames
using NLsolve 
using Plots
using Roots

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
    filepath4 = joinpath(@__DIR__, "USBenchmark.txt")
    USBenchmark = CSV.read(filepath4, DataFrame)
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
    #ChinaBenchmark = DataFrame(load(joinpath(@__DIR__, "ChinaBenchmark.jld2")))
    ChinaBenchmark = CSV.read("C:\\Users\\peter\\.julia\\dev\\dev_econ_replication\\replication_files\\urban_accounting_welfare_replication\\ReplicationFiles\\MatlabPro\\ChinaBenchmark.txt",
    header=false, DataFrame) 
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

Nbar = sum(1000 .*shocks[:,5])
Nlar = maximum(1000 .*shocks[:,5])

N = 1000 .*shocks[:,5]

c1 = zeros(Float64, 193)
c2 = zeros(Float64, 193)

for i in 1:length(shocks[:,1])

    c1[i] = ubar + (1 + psi)*log((1 + psi))-psi*log(psi)
    c2[i] = (1-theta)*((exp(shocks[i,2])^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))

    shocks[i,3] = -log(c2[i])+((N[i]*((kappa^2)/pi))^(1/2))*(((1+psi)/c2[i])+(2/3)*exp(shocks[i,4])) + c1[i]
    shocks[i,3] = log(shocks[i,3])
end 

NbarNew = sum(N)

#Counterfactuals without Differences in One City Characteristic, κ = $(kappa), ω = $(eps), ζ = $(ksi)

# Counterfactual efficiency
count = 0

step = 0.1 + China

NbarNewTFP = 0

ubarTFP = ubar

NTFP = zeros(Float64, 193)

RTFP = 0

X0 = 0

while abs(Nbar-NbarNewTFP) > 0.02

    for i = 1:length(shocks[:,1])

        if China == 0
        
            global RTFP = sum(exp.(shocks[:,2])./((shocks[:,5].*1000).^eps).*(shocks[:,5].*1000 ./Nbar))

        else
            
            global RTFP = quantile(exp.(shocks[:,2])./((shocks[:,5].*1000).^eps), 0.49)

        end    

        
        if NbarNewTFP == 0
           global X0 = 1000*shocks[i,5]
        else
           global X0 = max(NTFP[i], 1)
        end

        a1 = ubarTFP + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(X0^ksi))

        a2 = (1-theta)*(((RTFP*(X0^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
        f(X) = (pi/(kappa^2))*((log(a2)-a1)/(((1+psi)/a2)+(2/3)*exp(shocks[i,4])))^2-X

        root = find_zero(f, X0)
    
        NTFP[i] = max(0, root)
        
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

println("NbarNewTFP = $NbarNewTFP") 
println("Maximum NTFP = $(maximum(NTFP))")
println("ubarTFP = $ubarTFP")



# Counterfactual ability
count = 0

step = 0.1 + China

NbarNewA = 0

ubarA = ubar

AA = 0

NA = zeros(Float64, 193)

X0 = 0

ATFP = 0

while abs(Nbar-NbarNewA) > 0.02

    for i = 1:length(shocks[:,1])
        
        if China == 0
        
            global AA = sum(exp.(shocks[:,3])./((shocks[:,5].*1000).^ksi).*shocks[:,5].*1000 ./Nbar)
            
        else
            
            global AA = quantile(exp.(shocks[:,3])./((shocks[:,5].*1000).^ksi), 0.49)
            
        end    
        
        if NbarNewA == 0
            global X0 = 1000*shocks[i,5]
        else
            global X0 = max(NA[i], 1)
        end
        
        global ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)  

        b1 = ubarA + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(X0^ksi)

        b2 = (1-theta)*(((ATFP*(X0^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
            
        g(X) = (pi/(kappa^2))*((log(b2)-b1)/(((1+psi)/b2)+(2/3)*exp(shocks[i,4])))^2-X

        root = find_zero(g, X0)

        NA[i] = max(0, root)
            
    end 
 
    global NbarNewA = sum(NA)
    global count = count + 1
    if count > 10
        global step = step/10
        global count = 0
    end
    
    if Nbar-NbarNewA > 0.02
        global ubarA = ubarA - step
    else
        global ubarA = ubarA + step
    end
    println(Nbar-NbarNewA)
end

println("NbarNewA = $NbarNewA") 
println("Maximum NA = $(maximum(NA))")
println("ubarA = $ubarA")


# Counterfactual EF
count = 0

step = 0.1

NbarNewEF = 0

ubarEF = ubar

NEF = zeros(Float64, 193)

AEF = 0

X0 = 0

ATFP = 0

while abs(Nbar-NbarNewEF) > 0.02

    for i = 1:length(shocks[:,1])
        
        if China == 0
        
            global AEF = sum(exp.(shocks[:,4]).*shocks[:,5].*1000 ./Nbar)
            
        else
            
            global AEF = quantile(exp.(shocks[:,4]), 0.49)
            
        end    
        
        if NbarNewEF == 0
            global X0 = 1000*shocks[i,5]
        else
            global X0 = max(NEF[i], 1)
        end
        
        global ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)

        d1 = ubarEF + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(X0^ksi))

        d2 = (1-theta)*(((ATFP*(X0^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
        z(X) = (pi/(kappa^2))*((log(d2)-d1)/(((1+psi)/d2)+(2/3)*AEF))^2-X

        root = find_zero(z, X0)
    
        NEF[i] = max(0, root)
        
    end 
 
    global NbarNewEF = sum(NEF)
    global count = count + 1
    if count > 10
        global step = step/10
        global count = 0
    end
    
    if Nbar-NbarNewEF > 0.02
        global ubarEF = ubarEF - step
    else
        global ubarEF = ubarEF + step
    end
    println(Nbar-NbarNewEF)
end

println("NbarNewEF = $NbarNewEF") 
println("Maximum NEF = $(maximum(NEF))") 
println("ubar EF = $ubarEF")


TMA = sum(max.((NA .- N), 0))/Nbar
TMTFP = sum(max.((NTFP .- N), 0))/Nbar
TMEF = sum(max.((NEF .- N), 0))/Nbar

SNR = sort(log.(1000 .*shocks[:,1]))
SN = sort(log.(N))
SNTFP = sort(log.(NTFP))
SNA = sort(log.(NA))
SNEF = sort(log.(NEF))

lprob = zeros(Float64, length(N))

for i = 1:length(N)
    lprob[i] = log((length(N) - i + 1)/length(N))
end


ubarTFP = trunc(ubarTFP, digits=4)
TMTFP = trunc(TMTFP, digits=4)
ubarA = trunc(ubarA, digits=4)
TMA = trunc(TMA, digits=4)
ubarEF = trunc(ubarEF, digits=4)
TMEF = trunc(TMEF, digits=4)

f1p1 = plot(SNR, lprob, color=:blue, linewidth=2, label="Actual", legend=:bottomleft, title="Model Utility = $(ubar)",
titlefontsize=8, xlabel="ln(population)", ylabel="ln(prob > population)", titlelocation=:left,
xlabelfontsize=8, ylabelfontsize=8, xlims=(11, 17))
plot!(f1p1, SN, lprob, label="Modeled", color=:red, legend=:bottomleft, xlims=(11, 17), 
linewidth=2)

f1p2 = plot(SN, lprob, label="Actual", legend=:bottomleft, 
title="Counterfactual Utility = $(ubarTFP), \n reallocation = $(TMTFP)", titlefontsize=8, titlelocation=:left,
xlabel="ln(population)",
ylabel="ln(prob > population)",
color=:blue, linewidth=2,
xlabelfontsize=8, ylabelfontsize=8, xlims=(11, 18))
plot!(f1p2, SNTFP, lprob, label="Avg. Efficiency", legend=:bottomleft, color=:green, linewidth=2, xlims=(11, 18))

f1p3 = plot(SN, lprob, label="Actual", legend=:bottomleft,
title="Counterfactual Utility = $(ubarA), \n reallocation = $(TMA)", titlefontsize=8, titlelocation=:left,
xlabel="ln(population)", ylabel="ln(prob > population)",
color=:blue, linewidth=2,
xlabelfontsize=8, ylabelfontsize=8, xlims=(11, 18))
plot!(f1p3, SNA, lprob, label="Avg. Amenities", legend=:bottomleft, color=:magenta, linewidth=2, xlims=(11, 18))

f1p4 = plot(SN, lprob, 
label="Actual", legend=:bottomleft,
xlabel="ln(population)",
ylabel="ln(prob > population)",
title="Counterfactual Utility = $ubarEF, \n reallocation = $(TMEF)", titlefontsize=8, titlelocation=:left,
color=:blue, linewidth=2,
xlabelfontsize=8, ylabelfontsize=8, xlims=(11, 17))
plot!(f1p4, SNEF, lprob, label="Avg. Exc. Frictions", legend=:bottomleft, color=:black, linewidth=2, xlims=(11, 17))

f1 = plot(f1p1, f1p2, f1p3, f1p4, layout=(2,2))
savefig(f1, "figure_1.png")

chA = ((NA .- N)./N)
chTFP = ((NTFP .- N)./N)
chEF = ((NEF .- N)./N)
PercChange = DataFrame(hcat(chA, chTFP, chEF), :auto)
ShocksModel = shocks[:, 1:4]

if China == 0 
    CSV.write("PercChangeUS.csv", PercChange)
    CSV.write("ShocksModelUS.csv", ShocksModel)
else
    save("PercChangeChina.csv", PercChange)
    save("ShocksModelChina.csv", ShocksModel)
end    
    

#Counterfactuals with Differences in Only One City Characteristic, κ = $(kappa), ω = $(eps), ζ = $(ksi)

#Only US
if China == 0
 
    # Counterfactual TFP only
    count = 0

    step = 0.1

    NbarNewTFPO = 0

    ubarTFPO = ubar

    NTFPO = zeros(Float64, 193)

    X0 = 0

    ATFP = 0

    while abs(Nbar-NbarNewTFPO) > 0.02

        for i = 1:length(shocks[:,1])
        
            if NbarNewTFPO == 0
                global X0 = 1000*shocks[i,5]
            else
                global X0 = max(NTFPO[i], 1)
            end
        
            global ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)

            e1 = ubarTFPO + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(X0^ksi)

            e2 = (1-theta)*(((ATFP*(X0^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
            v(X) = (pi/(kappa^2))*((log(e2)-e1)/(((1+psi)/(e2))+(2/3)*AEF))^2-X

            root = find_zero(v, X0)
    
            NTFPO[i] = max(0, root)
        
        end 
 
        global NbarNewTFPO = sum(NTFPO)
        global count = count + 1
        if count > 10
            global step = step/10
            global count = 0
        end
    
        if Nbar-NbarNewTFPO > 0.02
            global ubarTFPO = ubarTFPO - step
        else
            global ubarTFPO = ubarTFPO + step
        end
        println(Nbar-NbarNewTFPO)
    end

    println("NbarNewTFPO = $NbarNewTFPO") 
    println("Maximum NTFPO = $(maximum(NTFPO))") 
    println("ubarTFPO = $ubarTFPO")

    TMTFPO = sum(max.((NTFPO .- N), 0))/Nbar

    ubarTFPO = trunc(ubarTFPO, digits=4)
    TMTFPO = trunc(TMTFPO, digits=4)

    f2p1 = plot(SNR, lprob, label="Actual", legend=:bottomleft,
    xlabel="ln(population)", ylabel="ln(prob > population)",
    title="Model Utility = $(ubar)", titlefontsize=8, titlelocation=:left, 
    color=:red, linewidth=2,
    xlabelfontsize=8, ylabelfontsize=8, xlims=(11, 18))
    plot!(f2p1, SN, lprob, label="Modeled", legend=:bottomleft, color=:blue, linewidth=2)

    SNTFPO = sort(log.(NTFPO))

    f2p2 = plot(SN, lprob,
    xlabel="ln(population)", ylabel="ln(prob > population)", xlims=(11, 18),
    label="Actual",
    legend=:bottomleft,
    title="Counterfactual Utility = $(ubarTFPO), \n  Reallocation = $(TMTFPO)", titlefontsize=8, titlelocation=:left,
    color=:blue, linewidth=2,
    xlabelfontsize=8, ylabelfontsize=8)
    plot!(f2p2, SNTFPO, lprob, label="Efficiency Only", legend=:bottomleft, color=:green, linewidth=2)

    # Counterfactual Amenities Only
    count = 0

    step = 0.1

    NbarNewAO = 0

    ubarAO = ubar

    NAO = zeros(Float64, 193)

    X0 = 0

    while abs(Nbar-NbarNewAO) > 0.02

        for i = 1:length(shocks[:,1])
        
            if NbarNewAO == 0
                global X0 = 1000*shocks[i,5]
            else
                global X0 = max(NAO[i], 1)
            end
        
            g1 = ubarAO + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(X0^ksi))

            g2 = (1-theta)*(((RTFP*(X0^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
                
            s(X) = (pi/(kappa^2))*((log(g2)-g1)/(((1+psi)/g2)+(2/3)*AEF))^2-X

            root = find_zero(s, X0)
    
            NAO[i] = max(0, root)
        
        end 
 
        global NbarNewAO = sum(NAO)
        global count = count + 1
        if count > 10
            global step = step/10
            global count = 0
        end
    
        if Nbar-NbarNewAO > 0.02
            global ubarAO = ubarAO - step
        else
            global ubarAO = ubarAO + step
        end
        println(Nbar-NbarNewAO)
    end

    println("NbarNewAO = $NbarNewAO") 
    println("Maximum NAO = $(maximum(NAO))") 
    println("ubarAO = $ubarAO")

    TMAO = sum(max.((NAO .- N), 0))/Nbar

    SNAO = sort(log.(NAO))

    ubarAO = trunc(ubarAO, digits=4)
    TMAO = trunc(TMAO, digits=4)

    f2p3 = plot(SN, lprob, 
    xlabel="ln(population)", ylabel="ln(prob > population)", xlims=(11, 18),
    legend=:bottomleft,
    label="Actual",
    title="Counterfactual Utility = $(ubarAO), \n Reallocation = $(TMAO)", titlefontsize=8, titlelocation=:left,
    color=:blue, linewidth=2,
    xlabelfontsize=8, ylabelfontsize=8)
    plot!(f2p3, SNAO, lprob, label="Amenities Only", legend=:bottomleft, color=:magenta, linewidth=2)

    # Counterfactual Excessive Frictions Only
    count = 0

    step = 0.1

    NbarNewEFO = 0

    ubarEFO = ubar

    NEFO = zeros(Float64, 193)

    X0 = 0

    while abs(Nbar-NbarNewEFO) > 0.02

        for i = 1:length(shocks[:,1])
        
            if NbarNewEFO == 0
                global X0 = 1000*shocks[i,5]
            else
                global X0 = max(NEFO[i], 1)
            end
        
            h1 = ubarEFO + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(X0^ksi)

            h2 = (1-theta)*(((RTFP*(X0^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
            
            k(X) = (pi/(kappa^2))*((log(h2)-h1)/(((1+psi)/h2)+(2/3)*exp(shocks[i,4])))^2-X

            root = find_zero(k, X0)
    
            NEFO[i] = max(0, root)
        
        end 
 
        global NbarNewEFO = sum(NEFO)
        global count = count + 1
        if count > 10
            global step = step/10
            global count = 0
        end
    
        if Nbar-NbarNewEFO > 0.02
            global ubarEFO = ubarEFO - step
        else
            global ubarEFO = ubarEFO + step
        end
        println(Nbar-NbarNewEFO)
    end

    println("NbarNewEFO  = $NbarNewEFO") 
    println("Maximum NEFO = $(maximum(NEFO))") 
    println("ubarEFO = $ubarEFO")

    TMEFO = sum(max.((NEFO .- N), 0))/Nbar

    SNEFO = sort(log.(NEFO))

    ubarEFO = trunc(ubarEFO, digits=4)
    TMEFO = trunc(TMEFO, digits=4)

    f2p4 = plot(SN, lprob, 
    label="Actual", xlabel="ln(population)", ylabel="ln(prob > population)", legend=:bottomleft,
    title="Counterfactual Utility = $(ubarEFO), \n Reallocation = $(TMEFO)", titlefontsize=8, titlelocation=:left,
    color=:blue, linewidth=2,
    xlabelfontsize=8, ylabelfontsize=8)
    plot!(f2p4, SNEFO, lprob, label="Exc. Frictions Only", legend=:bottomleft, color=:black, linewidth=2)

    f2 = plot(f2p1, f2p2, f2p3, f2p4, layout=(2,2))
    savefig(f2, "figure_2.png")

    # Counterfactual Without Shocks
    count = 0

    step = 0.1

    NbarNewNS = 0

    ubarNS = ubar

    NNS = zeros(Float64, 193)

    X0 = 0

    while abs(Nbar-NbarNewNS) > 0.02

        for i = 1:length(shocks[:, 1])
        
            if NbarNewNS == 0
                global X0 = 1000*shocks[i,5]
            else
                global X0 = max(NNS[i], 1)
            end
        
            j1 = ubarNS + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(X0^ksi)

            j2 = (1-theta)*(((RTFP*(X0^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
            
            q(X) = (pi/(kappa^2))*((log(j2)-j1)/(((1+psi)/j2)+(2/3)*AEF))^2-X

            root = find_zero(q, X0)
    
            NNS[i] = max(0, root)
        
        end 
 
        global NbarNewNS = sum(NNS)
        global count = count + 1
        if count > 10
            global step = step/10
            global count = 0
        end
    
        if Nbar-NbarNewNS > 0.02
            global ubarNS = ubarNS - step
        else
            global ubarNS = ubarNS + step
        end
        println(Nbar-NbarNewNS)
    end

    println("NbarNewNS = $NbarNewNS") 
    println("Maximum NNS = $(maximum(NNS))") 
    println("ubarNS = $ubarNS")

    SNNS = sort(log.(NNS))

    ubarNS = trunc(ubarNS, digits=4)

    f3 = plot(SN, lprob, 
    label="Actual", legend=:bottomleft,
    title="Counterfactual Utility = $(ubarNS)", titlefontsize=8, xlims=(11, 18), titlelocation=:left,
    xlabel="ln(population)", ylabel="ln(prob > population)",
    color=:blue, linewidth=2,
    xlabelfontsize=8, ylabelfontsize=8)
    plot!(f3, SNNS, lprob, label="No Shocks", legend=:bottomleft, color=:black, linewidth=2)
    savefig(f3, "figure_3.png")
end 