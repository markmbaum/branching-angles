using DrWatson
@quickactivate "Martian Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Plots
using Shapefile
using Serialization
using DataFrames
using UnPack

##

assembly = deserialize(datadir("exp_pro", "conus_assembly_serialized"));
angles = deserialize(datadir("exp_pro", "conus_angles_serialized"));

@unpack geoms, orders, networks, junctions, neighbors = assembly;

##

function plotjpt!(p, junction)::Nothing
    plot!(p, [junction[1]], [junction[2]],
        markershapes=:circle,
        markercolors=:black,
        markersize=4,
        legend=false)
    return nothing
end

function plotvalley!(p, geom, order, maxorder, junction)::Nothing
    #valley line
    x = [pt.x for pt ∈ geom.points]
    y = [pt.y for pt ∈ geom.points]
    plot!(p, x, y,
        color=:blue,
        alpha=1/(maxorder - order + 1),
        linewidth=order^0.85)
    #stream order
    m = length(x) ÷ 2
    annotate!(p, x[m], y[m], ("$order", 9, :red, :center))
    #angle regression
    θ = valleyangle(x, y, junction...)
    s = tan(θ)
    Δx = maximum(x) - minimum(x)
    Δy = maximum(y) - minimum(y)
    if abs(s*Δx) > Δy
        Δx *= Δy/abs(s*Δx)
    end
    x₁, y₁ = junction
    if -π/2 <= θ <= π/2
        x₂ = x₁ + Δx
        y₂ = y₁ + s*Δx
    else
        x₂ = x₁ - Δx
        y₂ = y₁ - s*Δx
    end
    plot!(p, [x₁, x₂], [y₁, y₂],
        color=:black,
        linestyle=:dash)
    return nothing
end

function plotjunction(idx::Int, geoms, orders, junctions, neighbors)
    J = junctions[idx]
    I = neighbors[idx]
    G = geoms[I]
    O = orders[I]
    p = plot(; aspect_ratio=:equal)
    plotjpt!(p, J)
    for i ∈ 1:length(G)
        println(i)
        plotvalley!(p, G[i], O[i], maximum(O), J)
    end
    B = findbranchingangles(G, J, O, I)
    println(B)
    null = isempty(B) ? " (empty)" : ""
    title!(p, "Case Number: $(B.case)$null")
    return p
end

function plotjunction(geoms, orders, junctions, neighbors)
    idx = rand(1:length(junctions))
    plotjunction(idx, geoms, orders, junctions, neighbors)
end

##

case2 = findall(b -> b.case == 2, angles)

##

plotjunction(geoms, orders, junctions, neighbors)