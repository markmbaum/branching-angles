using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Arrow
using CSV
using DataFrames

##

fn = datadir("exp_pro", "conus_angles_raw.csv")

##

df = fn |> CSV.File |> DataFrame

##

for col ∈ names(df)
    t = eltype(df[!, col])
    if t == Float64
        df[!, col] = Float32.(df[!, col])
    elseif t == Int64
        df[!, col] = Int32.(df[!, col])
    end
end

##

Arrow.write(replace(fn, ".csv" => ".feather"), df)
