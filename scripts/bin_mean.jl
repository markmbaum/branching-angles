using DrWatson
@quickactivate "Martian Branching Angles"
using DataFrames
using CSV
using Statistics: mean

##

function binmean(x::AbstractVector, y::AbstractVector; bins::Int=8, null::Function=isnan)
    i = Not(null.(x))
    x = x[i]
    y = y[i]
    b = collect(LinRange(minimum(x), maximum(x), bins+1))
    println(b)
    c = (b[1:end-1] .+ b[2:end])/2
    m = [mean(y[(x .>= b[i]) .& (x .<= b[i+1])]) for i ∈ 1:bins]
    return c, m
end

##

df = CSV.read(datadir("exp_pro", "conus_angles.csv"), DataFrame)