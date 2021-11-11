using DrWatson
@quickactivate "Martian Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Serialization
using DataFrames
using CSV
using Shapefile
using LinearAlgebra: ⋅ #dot product
using MultiAssign
using Base.Threads: @threads

#mean radius of Mars
const 𝐑ₘ = 3.398e6

##

function latlon2sph(lat, lon)
    θ = (-lat/180)*π + π/2
    ϕ = (lon/180)*π + π
    return θ, ϕ
end

##

#branching angles are already calculated
angles = deserialize(datadir("dir_pro", "mars_branching_angles"))
#create a dataframe from the angles
df = DataFrame(angles)
#use the unprojected (geographic coords) file for subsequent slope calcs
fn = datadir("exp_raw", "mars-networks", "Hynek_Valleys_geology_draped.shp")
db = DataFrame(Shapefile.Table(fn))

## fill in slope for each valley

L = size(db)[1]
db[!,:slope] = zeros(L)
@threads for i = 1:L
    geom = db[i,:geometry]
    P = geom.points
    N = length(P)
    lon = [p.x for p ∈ P]
    lat = [p.y for p ∈ P]
    #center on the mean
    lon .-= sum(lon)/N
    lat .-= sum(lat)/N
    #convert to true spherical coordinates
    @multiassign θ, ϕ = zeros(N)
    for i = 1:N
        θ[i], ϕ[i] = latlon2sph(lat[i], lon[i])
    end
    #equirectangular projection
    x = 𝐑ₘ*ϕ
    y = 𝐑ₘ*(π .- θ)
    #elevation is in the "measure" vector
    z = geom.measures
    #orthogonal distance regression for slope
    db[i,:slope] = valleyslope(x, y, z)
end

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

L = size(df)[1]
@multiassign df[!,:lat], df[!,:lon] = zeros(L)
@threads for i = 1:L
    idx₁ = df[i,"index A"]
    idx₂ = df[i,"index B"]
    geom₁ = db[idx₁,:geometry]
    geom₂ = db[idx₂,:geometry]
    df[i,:lon], df[i,:lat] = intersection(geom₁.points, geom₂.points)
end

## write the finalized dataframe to csv
CSV.write(datadir("exp_pro", "mars_branching_angles.csv"), df)