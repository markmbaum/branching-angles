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

##

#mean radius of Mars
const 𝐑ₘ = 3.398e6

function projgeom(geom)
    P = geom.points
    lon = [p.x for p ∈ P]
    lat = [p.y for p ∈ P]
    #recenter
    lon .-= (maximum(lon) + minimum(lon))/2
    lat .-= (maximum(lat) + minimum(lat))/2
    #convert to radians
    θ = lat*(π/180)
    ϕ = lon*(π/180)
    #convert to approximate distance coordinates [m]
    x = 𝐑ₘ*ϕ
    y = 𝐑ₘ*θ
    return x, y
end

##

#branching angles are already calculated in "find_branching_angles.jl" script
df = deserialize(
    datadir(
        "exp_pro",
        "mars_angles_serialized"
    )
) |> DataFrame
#use the unprojected (geographic coords) file for subsequent slope calcs
db = datadir(
    "exp_raw",
    "mars-networks",
    "Hynek_Valleys_geology_draped.shp"
) |> Shapefile.Table |> DataFrame

## fill in slope for each valley

L = size(db, 1)
db[!,:slope] = fill(NaN, L)
@threads for i = 1:10#L
    geom = db[i,:geometry]
    #project into appropriate distance coordinates
    x, y = projgeom(geom)
    #elevation is in the "measure" vector
    z = geom.measures
    #orthogonal distance regression for slope
    db[i,:slope] = valleyslope(x, y, z)
end

## DROP CASES 1,2,3

df = dropcases(df, [1,2,3])

##

#indices of streams for all angles
A = df[:,"index A"]
B = df[:,"index B"]

#map attributes into the angles dataframe
df[!,"geology A"] = db[A,:Unit]
df[!,"geology B"] = db[B,:Unit]
df[!,"slope A"] = db[A,:slope]
df[!,"slope B"] = db[B,:slope]

## fill in junction locations in latitude and longitude

L = size(df, 1)
@multiassign df[!,:lat], df[!,:lon] = fill(NaN, L)
@threads for i = 1:L
    idx₁ = df[i,"index A"]
    idx₂ = df[i,"index B"]
    geom₁ = db[idx₁,:geometry]
    geom₂ = db[idx₂,:geometry]
    df[i,:lon], df[i,:lat] = intersection(geom₁.points, geom₂.points)
end

## write the finalized dataframe to csv
CSV.write(datadir("exp_pro", "mars_branching_angles.csv"), df)