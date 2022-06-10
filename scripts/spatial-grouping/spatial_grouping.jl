using DrWatson
@quickactivate "Branching Angles"

using Pkg
Pkg.instantiate()

push!(LOAD_PATH, srcdir())
using BranchingAngles
using Base.Threads: @threads, nthreads
using Arrow
using DataFrames
using StatsBase: zscore, mean, std

## load and prepare CONUS data

df = datadir("exp_pro", "conus_angles.feather") |> Arrow.Table |> DataFrame
df = mapcols(Vector, df)
#computes logslope, minorder, maxorder
derivedcols!(df)
#renames precip and temperature columns
renamePT!(df)
#column selection
cols = [:angle, :P, :T, :AI, :EVI, :SSM, :logslope, :minorder, :maxorder]
#selection & standardization
select!(df, vcat(cols, :x, :y))
transform!(df, cols .=> zscore .=> cols)

## function for spatial grouping and averaging

function neighborhood(df, ngroups::Int, Δ, minsamples::Int)::DataFrame
    x = df.x .|> Float32
    y = df.y .|> Float32
    cols = filter(c -> (c != "x") & (c != "y"), names(df))
    M = df[!,cols] |> Matrix .|> Float32
    b = zeros(Bool, size(M,1))
    ne = DataFrame(zeros(Float32, ngroups, length(cols)), cols)
    ne[!,:n] = zeros(Int32, ngroups)
    ne[!,:Δ] = fill(Δ, ngroups)
    g = 0
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    while g < ngroups
        #pick a random spot
        xᵣ = (xmax - xmin)*rand() + xmin
        yᵣ = (ymax - ymin)*rand() + ymin
        #clear the boolean vector
        b .= false
        #find observations inside the box
        for i ∈ 1:length(b)
            if (xᵣ - Δ/2 < x[i] < xᵣ + Δ/2)
                if (yᵣ - Δ/2 < y[i] < yᵣ + Δ/2)
                    b[i] = true
                end
            end
        end
        #check if there are enough observations in the 'hood
        samples = sum(b)
        if samples >= minsamples
            g += 1
            ne[g,:n] = samples
            m = mapcols(mean, df[b,:])
            for c ∈ cols
                ne[g,c] = m[1,c]
            end
        end
    end
    return ne
end

##

#define box sizes
xmin, xmax = extrema(df.x)
Δ = LinRange((xmax - xmin)/50, (xmax - xmin)/10, nthreads())

#compute in parallel
nes = Vector{DataFrame}(undef, nthreads())
@threads for i ∈ 1:nthreads()
    nes[i] = neighborhood(df, 100, Δ[i], 100)
end

#concatenate and save
Arrow.write(
    datadir(
        "exp_pro",
        "spatial_grouping.feather"
    ),
    vcat(nes...)
)