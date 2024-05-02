using CSV 
using DataFrames
using NLsolve 
using Plots
using Roots
using Statistics

#Select Country 
China = 1 #Make 0 if US and 1 if China

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
    filepath45 = joinpath(@__DIR__, "ChinaBenchmark.txt")
    ChinaBenchmark = CSV.read(filepath45, DataFrame) 
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

c1 = zeros(Float64, 212)
c2 = zeros(Float64, 212)

for i in 1:length(shocks[:,1])

    c1[i] = ubar + (1 + psi)*log((1 + psi))-psi*log(psi)
    c2[i] = (1-theta)*((exp(shocks[i,2])^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))

    shocks[i,3] = -log(c2[i])+((N[i]*((kappa^2)/pi))^(1/2))*(((1+psi)/c2[i])+(2/3)*exp(shocks[i,4])) + c1[i]
    shocks[i,3] = log(shocks[i,3])
end 

NbarNew = sum(N)

# Counterfactual EF 50th Percentile
count = 0

step = 0.1

NbarNewEF = 0

ubarEF = ubar

NEF = zeros(Float64, 212)

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


TMEF = sum(max.((NEF .- N), 0))/Nbar
SN = sort(log.(N))
SNEF = sort(log.(NEF))
lprob = zeros(Float64, length(N))

for i = 1:length(N)
    lprob[i] = log((length(N) - i + 1)/length(N))
end

#Counterfactual EF 10th Percentile 

count = 0

step = 0.1

NbarNewEF = 0

ubarEF3 = ubar

NEF = zeros(Float64, 212)

AEF = 0

X0 = 0

ATFP = 0

while abs(Nbar-NbarNewEF) > 0.02

    for i = 1:length(shocks[:,1])
        
        if China == 0
        
            global AEF = sum(exp.(shocks[:,4]).*shocks[:,5].*1000 ./Nbar)
            
        else
            
            global AEF = quantile(exp.(shocks[:,4]), 0.09)
            
        end    
        
        if NbarNewEF == 0
            global X0 = 1000*shocks[i,5]
        else
            global X0 = max(NEF[i], 1)
        end
        
        global ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)

        d1 = ubarEF3 + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(X0^ksi))

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
        global ubarEF3 = ubarEF3 - step
    else
        global ubarEF3 = ubarEF3 + step
    end
    println(Nbar-NbarNewEF)
end

println("NbarNewEF = $NbarNewEF") 
println("Maximum NEF = $(maximum(NEF))") 
println("ubar EF = $ubarEF3")


TMEF = sum(max.((NEF .- N), 0))/Nbar
SN = sort(log.(N))
SNEF3 = sort(log.(NEF))
lprob = zeros(Float64, length(N))

for i = 1:length(N)
    lprob[i] = log((length(N) - i + 1)/length(N))
end

#Counterfactual EF 90th Percentile 

count = 0

step = 0.1

NbarNewEF = 0

ubarEF2 = ubar

NEF = zeros(Float64, 212)

AEF = 0

X0 = 0

ATFP = 0

while abs(Nbar-NbarNewEF) > 0.02

    for i = 1:length(shocks[:,1])
        
        if China == 0
        
            global AEF = sum(exp.(shocks[:,4]).*shocks[:,5].*1000 ./Nbar)
            
        else
            
            global AEF = quantile(exp.(shocks[:,4]), 0.89)
            
        end    
        
        if NbarNewEF == 0
            global X0 = 1000*shocks[i,5]
        else
            global X0 = max(NEF[i], 1)
        end
        
        global ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)

        d1 = ubarEF2 + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(X0^ksi))

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
        global ubarEF2 = ubarEF2 - step
    else
        global ubarEF2 = ubarEF2 + step
    end
    println(Nbar-NbarNewEF)
end

println("NbarNewEF = $NbarNewEF") 
println("Maximum NEF = $(maximum(NEF))") 
println("ubar EF = $ubarEF2")


TMEF = sum(max.((NEF .- N), 0))/Nbar
SN = sort(log.(N))
SNEF2 = sort(log.(NEF))
lprob = zeros(Float64, length(N))

for i = 1:length(N)
    lprob[i] = log((length(N) - i + 1)/length(N))
end

#Changing Excessive Frictions in China 

ubarEF = trunc(ubarEF, digits=4)
ubarEF2 = trunc(ubarEF2, digits=4)
ubarEF3 = trunc(ubarEF3, digits=4)

f1 = plot(SNEF, lprob, xlabel="ln(population)", ylabel="ln(prob > population)", color=:black, linewidth=2,
xlabelfontsize=8, ylabelfontsize=8, xlims=(11, 17), label="All excessive frictions at 50th percentile, utility=$(ubarEF)", 
legend=:bottomleft)
plot!(f1, SNEF2, lprob, color=:red, linewidth=2, label="All excessive frictions at 90th percentile, utility=$(ubarEF2)")
plot!(f1, SNEF3, lprob, color=:blue, linewidth=2, label="All excessive frictions at 10th percentile, utility=$(ubarEF3)")

savefig(f1, "figure_9.png")