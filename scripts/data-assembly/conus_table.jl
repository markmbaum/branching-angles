using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
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
df = CSV.read(datadir("exp_pro", "conus_angles_raw.csv"), DataFrame)
#big file, takes some time to load
db =
    datadir("exp_raw", "conus-networks", "conus_networks_lambert.shp") |>
    Shapefile.Table |>
    DataFrame

## DROP CASES 1,2,3

df = dropcases(df, 1, 2, 3)

## drop the non-stream/river identifiers

#shape types
FTYPEA = db[df[:, "index A"], :FTYPE]
FTYPEB = db[df[:, "index B"], :FTYPE]

df = df[(FTYPEA.=="StreamRiver").&(FTYPEB.=="StreamRiver"), :]

## add slope columns

df[!, "slope A"] = conditionslope.(db[df[!, "index A"], "SLOPE"])
df[!, "slope B"] = conditionslope.(db[df[!, "index B"], "SLOPE"])

## write the finalized dataframe to csvs

CSV.write(datadir("exp_pro", "conus_anglesX.csv"), df)
