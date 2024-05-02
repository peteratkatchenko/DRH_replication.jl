module DRH_replication

export drh_replicate

using StatFiles
using DataFrames
using Econometrics
using CSV
using JLD2
using Statistics
using GLM
using NLsolve
using Plots
using Roots

function drh_replicate()

    include(joinpath(@__DIR__, "InputsNumericalUS.jl"))

    include(joinpath(@__DIR__, "InputsNumericalChina.jl"))

    include(joinpath(@__DIR__, "UAWMainUSA23.jl"))

    include(joinpath(@__DIR__, "UAWMainUSA4.jl"))

    include(joinpath(@__DIR__, "UAWMainChina8.jl"))

    include(joinpath(@__DIR__, "UAWMainChina9.jl"))

end 


end #End of the main development replication file 