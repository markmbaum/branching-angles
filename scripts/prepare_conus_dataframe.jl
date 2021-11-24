using DrWatson
@quickactivate "Martian Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Serialization
using DataFrames
using CSV
using Shapefile
using MultiAssign

##

#handles issues with the slopes in the CONUS data
#  null values are -9998.0 for some reason
conditionslope(s::Real)::Float64 = (s == -9998.0) ? NaN : Float64(s)

##

#branching angles are already calculated in "find_branching_angles.jl" script
df = deserialize(
    datadir(
        "exp_pro",
        "conus_angles_serialized"
    )
) |> DataFrame
#big file, takes some time to load
db = datadir(
    "exp_raw",
    "conus-networks",
    "NHDFlowline_Network_NoMZ_w_basins.shp"
) |> Shapefile.Table |> DataFrame

## DROP CASES 1,2,3

df = dropcases(df, [1,2,3])

## drop the non-stream/river identifiers

#shape types
FTYPEA = db[df[:,"index A"],:FTYPE]
FTYPEB = db[df[:,"index B"],:FTYPE]

df = df[(FTYPEA .== "StreamRiver") .& (FTYPEB .== "StreamRiver"), :]

## add slope columns

@multiassign df[!,"slope A"], df[!,"slope B"] = fill(NaN, size(df,1))
for i ∈ 1:size(df,1)
    df[i,"slope A"] = conditionslope(db[df[i,"index A"], "SLOPE"])
    df[i,"slope B"] = conditionslope(db[df[i,"index B"], "SLOPE"])
end

## write the finalized dataframe to csv

CSV.write(datadir("exp_pro", "conus_branching_angles.csv"), df)