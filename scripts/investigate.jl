using DrWatson
@quickactivate "Martian Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Serialization
using Shapefile
using DataFrames
using PyPlot

pygui(true)

##

function plotgeom(geom)
    plot(
        [p.x for p ∈ geom.points],
        [p.y for p ∈ geom.points]
    )
end

function plotgeoms(geoms)
    for geom ∈ geoms
        plotgeom(geom)
    end
end

##

df = datadir(
    "exp_raw",
    "mars-networks",
    "Hynek_Valleys_geology.shp"
) |> Shapefile.Table |> DataFrame

##

fid = [8253, 8224, 8285, 8221, 8296, 8481, 8328, 8292, 7997, 8204, 8215, 7808, 7816, 7858]

idx = [findfirst(x->x==i,df[!,:FID_1]) for i ∈ sort(fid)]

NA = assemblenetworks(
    df[idx,:geometry],
    df[idx,:VNOrder]
)

println(NA.networks)
println(NA.neighbors)
for J ∈ NA.junctions
    println("$(J.x) $(J.y)")
end
plotgeoms(df[idx,:geometry])

branchingangles(NA) |> DataFrame
