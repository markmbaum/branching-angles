using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using DataFrames
using CSV
using Statistics
using Random

## read and rearrange CONUS table

df = CSV.File(
    datadir(
        "exp_pro",
        "conus_angles.csv"
    )
) |> DataFrame

df[!,:P] = df[!,:ppt_annual]
df[!,:T] = df[!,:tmean_annual]
df[!,:logslope] = log10.(df[!,"slope A"])/2 .+ log10.(df[!,"slope B"])
df[!,:maxorder] = maximum.(zip(df[!,"order A"], df[!,"order B"]))
df[!,:minorder] = minimum.(zip(df[!,"order A"], df[!,"order B"]))

cols = Symbol[
    :angle,
    :P,
    :AI,
    :NDVI,
    :EVI,
    :SSM,
    :SUSM,
    :SMP,
    :logslope,
    :maxorder,
    :minorder
]

df = df[:,[cols; :x; :y]]

##

bounds(X::AbstractVector) = (minimum(X), maximum(X))

function findbox(xᵢ, yᵢ, x, y, Δx, Δy)
    idx = Set{Int64}()
    @assert length(x) == length(y)
    @inbounds for j ∈ 1:length(x)
        if (abs(x[j] - xᵢ) < Δx/2) & (abs(y[j] - yᵢ) < Δy/2)
            push!(idx, j)
        end
    end
    return collect(idx)
end

function geosample(df, xcol, ycol, Δx, Δy, N, f::F=mean, seed=1) where {F}
    #new frame for neighborhood samples
    gs = DataFrame(cols .=> fill(zeros(N), length(cols)))
    gs[!,:n] = zeros(Int64, N)
    #randomly select
    x, y = df[!,xcol], df[!,ycol]
    xmin, xmax = bounds(x)
    ymin, ymax = bounds(y)
    rng = Xoshiro(seed)
    i = 1
    while i <= N
        xᵢ = (xmax - xmin)*rand(rng) + xmin
        yᵢ = (ymax - ymin)*rand(rng) + ymin
        sl = view(df, findbox(xᵢ, yᵢ, x, y, Δx, Δy), 1:size(df,2)-2)
        n = size(sl,1)
        if n > 10
            gs[i,1:end-1] = mapcols(f, sl)[1,:]
            gs[i,:n] = n
            i += 1
        end
    end
    return gs
end

println("x bounds = $(bounds(df.x))\ny bounds = $(bounds(df.y))")

gs = geosample(df, :x, :y, 1e5, 1e5, 10000)