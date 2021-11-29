using DrWatson
@quickactivate "Martian Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Serialization
using DataFrames
using CSV
using Shapefile
using MultiAssign
using Base.Threads: @threads

## functions

#handles issues with the slopes in the CONUS data
#  null values are -9998.0 for some reason
conditionslope(s::Real)::Float64 = (s == -9998.0) ? NaN : Float64(s)

##

#branching angles are already calculated in "find_branching_angles.jl" script
df = deserialize(
    datadir(
        "exp_pro",
        "conus_angles"
    )
) |> DataFrame
#big file, takes some time to load
db = datadir(
    "exp_raw",
    "conus-networks",
    "conus_networks_lambert.shp"
) |> Shapefile.Table |> DataFrame

## DROP CASES 1,2,3

df = dropcases(df, 1, 2, 3)

## drop the non-stream/river identifiers

#shape types
FTYPEA = db[df[:,"index A"],:FTYPE]
FTYPEB = db[df[:,"index B"],:FTYPE]

df = df[(FTYPEA .== "StreamRiver") .& (FTYPEB .== "StreamRiver"), :]

## indices for adding more columns

A = df[!,"index A"]
B = df[!,"index B"]

## add slope columns

df[!,"slope A"] = conditionslope.(db[A,"SLOPE"])
df[!,"slope B"] = conditionslope.(db[B,"SLOPE"])

## fill in junction locations in projected coords

@multiassign df[!,:x], df[!,:y] = fill(NaN, size(df, 1))
@threads for i = 1:size(df, 1)
    geom₁ = db[df[i,"index A"],:geometry]
    geom₂ = db[df[i,"index B"],:geometry]
    x, y = intersection(geom₁.points, geom₂.points)
    df[i,:x] = x
    df[i,:y] = y
end

## write the finalized dataframe to csv

CSV.write(datadir("exp_pro", "conus_angles_initial.csv"), df)