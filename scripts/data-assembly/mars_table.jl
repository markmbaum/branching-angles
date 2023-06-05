using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using DataFrames
using CSV
using Shapefile
using MultiAssign
using Base.Threads: @threads

##

#mean radius of Mars
const 𝐑ₘ = 3.398e6

#for slope calculations
function projectgeometry(geom)
    P = geom.points
    lon = [p.x for p ∈ P]
    lat = [p.y for p ∈ P]
    #recenter
    lon .-= (maximum(lon) + minimum(lon)) / 2
    lat .-= (maximum(lat) + minimum(lat)) / 2
    #convert to radians
    θ = lat * (π / 180)
    ϕ = lon * (π / 180)
    #convert to approximate distance coordinates [m]
    x = 𝐑ₘ * ϕ
    y = 𝐑ₘ * θ
    return x, y
end

##

#branching angles are already calculated in "branching_angles.jl" script
df = CSV.read(datadir("exp_pro", "mars_angles_raw.csv"), DataFrame)
#need the unprojected (geographic coords) file for subsequent slope calcs
db =
    datadir("exp_raw", "mars-networks", "Hynek_Valleys_geology_draped.shp") |>
    Shapefile.Table |>
    DataFrame

## DROP CASES 1,2,3

df = dropcases(df, 1, 2, 3)

## fill in slope for every valley

L = size(db, 1)
db[!, :slope] = fill(NaN, L)
@threads for i = 1:L
    geom = db[i, :geometry]
    #project into appropriate distance coordinates
    x, y = projectgeometry(geom)
    #elevation is in the M vector
    z = geom.measures
    #orthogonal distance regression for slope
    db[i, :slope] = valleyslope(x, y, z)
end

##

#map slope into the angles dataframe
df[!, "slope A"] = db[df[:, "index A"], :slope]
df[!, "slope B"] = db[df[:, "index B"], :slope]

## fill in junction locations in latitude and longitude

@multiassign df[!, :lat], df[!, :lon] = fill(NaN, size(df, 1))
@threads for i = 1:size(df, 1)
    geom₁ = db[df[i, "index A"], :geometry]
    geom₂ = db[df[i, "index B"], :geometry]
    p = endintersection(
        geom₁.points[1],
        geom₁.points[end],
        geom₂.points[1],
        geom₂.points[end],
    )
    df[i, :lon] = p.x
    df[i, :lat] = p.y
end

## write the finalized dataframe to csv

CSV.write(datadir("exp_pro", "mars_angles.csv"), df)
