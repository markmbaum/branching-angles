using DrWatson
@quickactivate "Martian Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using DataFrames
using CSV

##

function run(fnshp::String, orders::Symbol, fnout::String)::Nothing
    #network and branch into a dataframe of unvarnished results
    df = assemblenetworks(fnshp, orders) |> branchingangles |> DataFrame
    #write to a csv file
    CSV.write(fnout, df)
    return nothing
end

## Martian streams

run(
    datadir(
        "exp_raw",
        "mars-networks",
        "Hynek_Valleys_geology_draped_mercator.shp"
    ),
    :VNOrder,
    datadir(
        "exp_pro",
        "mars_angles_raw.csv"
    )
)

## contiguous United States (CONUS) streams

run(
    datadir(
        "exp_raw",
        "conus-networks",
        "conus_networks_lambert.shp"
    ),
    :StreamOrde, #[sic]
    datadir(
        "exp_pro",
        "conus_angles_raw.csv"
    )
)
