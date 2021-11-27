module BranchingAngles

using Base.Threads: @spawn, Task, fetch
using Shapefile
using Optim
using ProgressMeter: @showprogress
using PrettyTables
using UnPack
using DataFrames: DataFrame
import DataFrames
using MultiAssign

#------------------------------------------------------------------------------
export intersection, assemblenetworks
export NetworkAssembly

function intersectingends(p₁, p₂)::Bool
    #first check if the end points are intersecting (most likely by far)
    for a ∈ (p₁[1], p₁[end]), b ∈ (p₂[1], p₂[end])
        a == b && return true
    end
    return false
end

function intersection(p₁, p₂)
    for a ∈ (p₁[1], p₁[end]), b ∈ (p₂[1], p₂[end])
        a == b && return convert(Float64, a.x), convert(Float64, a.y)
    end
    error("no intersecting end points available")
    return NaN, NaN
end

function samepoint(a::NTuple{2,Float64}, b::Shapefile.Point)::Bool
    (a[1] == b.x) & (a[2] == b.y)
end

struct NetworkAssembly{T}
    geoms::Vector{T} #shapefile geometries
    orders::Vector{Int64} #stream order of each geometry
    networks::Vector{Vector{Int64}} #groups of stream/valley indices for each full network
    junctions::Vector{NTuple{2,Float64}} #all junction points
    neighbors::Vector{Vector{Int64}} #stream/valley indices meeting at each junction
    juncworks::Vector{Int64} #which network each junction belongs to
end

function Base.show(io::IO, A::NetworkAssembly{T}) where {T}
    println(io, "NetworkAssembly:")
    println(io, "  geometries of type $T")
    println(io, "  $(length(A.junctions)) junctions")
    println(io, "  $(length(A.networks)) networks")
    print(io, "  stream orders $(sort(unique(A.orders)))")
end

function assemblenetworks(fn::String, ordercol::Symbol)
    
    #load geometries from the target shapefile
    table = DataFrame(Shapefile.Table(fn))
    geoms = table[:,:geometry]
    orders = table[:,ordercol]
    
    #vector of vectors of indices, to group shapes
    networks = [[1]] #first shape automatically starts the first group
    #add every line to a group
    @showprogress 3000 "Assembling Networks " for i ∈ 2:length(geoms)
        #current shape
        geom = geoms[i]
        #intersection flag
        found = false
        #count backward because line is probably more likely to be in recent groups
        j = length(networks)
        while (j > 0) & !found
            k = length(networks[j])
            while (k > 0) & !found
                #geometry corresponding to group j, member k
                @inbounds geomⱼₖ = geoms[networks[j][k]]
                #test if the lines are intersecting
                if intersectingends(geom.points, geomⱼₖ.points)
                    #add geom to group j
                    push!(networks[j], i)
                    #break loops
                    found = true
                end
                k -= 1
            end
            j -= 1
        end
        #if ungrouped, start a new group
        !found && push!(networks, [i])
    end

    #find junction points and valleys meeting at each point
    junctions = NTuple{2,Float64}[]
    neighbors = Vector{Int64}[]
    juncworks = Int64[]
    @showprogress 1 "Finding junctions & neighbors " for i = 1:length(networks)
        network = networks[i]
        N = length(network)
        idx = length(junctions) #track current end of junction list
        if N > 1
            #find all junction points in the group
            for i ∈ 1:N
                geomᵢ = geoms[network[i]]
                for j ∈ i+1:N
                    geomⱼ = geoms[network[j]]
                    #see if there is an intersection point
                    if intersectingends(geomᵢ.points, geomⱼ.points)
                        #extract the intersection point
                        p = intersection(geomᵢ.points, geomⱼ.points)
                        #do not add a duplicate
                        if !(p ∈ junctions)
                            push!(junctions, p)
                            push!(juncworks, i)
                        end
                    end
                end
            end
            #now find which lines meet at the points
            for p ∈ @view junctions[idx+1:end]
                push!(neighbors, Int64[])
                for i ∈ network
                    p₁ = geoms[i].points[1]
                    p₂ = geoms[i].points[end]
                    if samepoint(p, p₁) | samepoint(p, p₂)
                        push!(neighbors[end], i)
                    end
                end
            end
        end
    end
    #final construction
    NetworkAssembly(geoms, orders, networks, junctions, neighbors, juncworks)
end

#------------------------------------------------------------------------------
export slope, valleyangle, valleyintersectionangle, valleyslope
export BranchingAngleResult
export branchinganglecases, branchingangles

distance(x₁, y₁, x₂, y₂) = sqrt((x₁ - x₂)^2 + (y₁ - y₂)^2)

orthogonaldistance(m, b, x, y) = abs(-m*x + y - b)/sqrt(m^2 + 1)

orthogonalloss(m, b, x, y, N::Int) = sum(orthogonaldistance(m, b, x[i], y[i])^2 for i ∈ 1:N)

struct ValleyCoordinates
    x::Vector{Float64}
    y::Vector{Float64}
    N::Int64
end

function ValleyCoordinates(x::AbstractArray, y::AbstractArray)
    @assert length(x) == length(y)
    N = length(x)
    #recenter the coordinates
    xm = (maximum(x) + minimum(x))/2
    ym = (maximum(y) + minimum(y))/2
    xc = collect(Float64, x)
    yc = collect(Float64, y)
    @inbounds for i = 1:N
        xc[i] -= xm
        yc[i] -= ym
    end
    #construct
    ValleyCoordinates(xc, yc, N)
end

#functor computes orthogonal distance loss function for fitting a line with slope m, intercept b
(O::ValleyCoordinates)(m, b) = orthogonalloss(m, b, O.x, O.y, O.N)
(O::ValleyCoordinates)(X) = O(X[1], X[2])

#use orthogonal distance to fit a line and return the slope
function ODRslope(x::AbstractArray, y::AbstractArray)::Float64
    #warehouse the coordinates with a convenient functor
    V = ValleyCoordinates(x, y)
    #intitial estimate of slope
    δx = V.x[end] - V.x[1]
    δy = V.y[end] - V.y[1]
    m₀ = δx == 0 ? 1e3 : δy/δx
    #fit for line parameters
    sol = optimize(
        V,
        [m₀, 0.0],
        Newton(),
        autodiff=:forward
    )
    return sol.minimizer[1]
end

function valleyslope(x::AbstractVector, y::AbstractVector, z::AbstractVector)::Float64
    @assert length(x) == length(y) == length(z)
    N = length(x)
    #compute cumulative horizontal distance
    d = distance.(view(x,1:N-1), view(y,1:N-1), view(x,2:N), view(y,2:N))
    c = zeros(N)
    cumsum!(view(c,2:N), d)
    #fit a line to get the slope estimate
    s = ODRslope(c, z)
    #positive values only
    return abs(s)
end

function valleyslope(geom)
    x = [p.x for p ∈ geom.points]
    y = [p.y for p ∈ geom.points]
    z = geom.measures
    valleyslope(x, y, z)
end

function valleyangle(x::AbstractVector, y::AbstractVector, xⱼ::Real, yⱼ::Real)
    @assert length(x) == length(y)
    N = length(x)
    #compute angle of coordinates
    θ = atan(ODRslope(x, y))
    #max distance between junction point and valley coordinates
    h = 0.0
    for i ∈ 1:N
        @inbounds d = distance(x[i], y[i], xⱼ, yⱼ)
        if d > h
            h = d
        end
    end
    #must decide which direction the valley is pointing
    p₁ = (xⱼ + h*cos(θ), yⱼ + h*sin(θ)) #test point 1
    p₂ = (xⱼ + h*cos(θ + π), yⱼ + h*sin(θ + π)) #test point 2
    d₁ = 0.0
    d₂ = 0.0
    @inbounds for i ∈ 1:N
        d₁ += distance(p₁[1], p₁[2], x[i], y[i])
        d₂ += distance(p₂[1], p₂[2], x[i], y[i])
    end
    #distance to test point 1 is larger, flip angle direction
    if d₁ > d₂
        θ = θ > 0 ? θ - π : θ + π
    end
    return θ
end

function valleyangle(geom, J)
    x = [p.x for p ∈ geom.points]
    y = [p.y for p ∈ geom.points]
    valleyangle(x, y, J...)
end

function minorangle(θ)
    while (θ > π) | (θ < -π)
        θ = (θ > π) ? θ - 2π : θ + 2π
    end
    return θ
end

function valleyintersectionangle(x₁::AbstractVector, y₁::AbstractVector, J₁::NTuple{2},
                                 x₂::AbstractVector, y₂::AbstractVector, J₂::NTuple{2})
    #direction angle of each valley
    θ₁ = valleyangle(x₁, y₁, J₁...)
    θ₂ = valleyangle(x₂, y₂, J₂...)
    #difference between the angles
    ϕ = θ₁ - θ₂
    #correct for results θ > π or θ < -π
    return abs(minorangle(ϕ))
end

function valleyintersectionangle(J, geom₁, geom₂)
    x₁ = [p.x for p ∈ geom₁.points]
    y₁ = [p.y for p ∈ geom₁.points]
    x₂ = [p.x for p ∈ geom₂.points]
    y₂ = [p.y for p ∈ geom₂.points]
    valleyintersectionangle(x₁, y₁, J, x₂, y₂, J)
end

#-----

struct BranchingAngleResult
    θ::Vector{Float64} #branching angles
    i::Vector{Int64} #indices of first shape in pairs
    j::Vector{Int64} #indices of second shape in pairs
    oᵢ::Vector{Int64} #stream orders of shapes i
    oⱼ::Vector{Int64} #stream orders of shapes j
    junction::NTuple{2,Float64} #junction point
    case::Int64 #branching case/type
    N::Int64 #vector lengths (number of angles)
end

Base.length(R::BranchingAngleResult) = R.N

function Base.show(io::IO, R::BranchingAngleResult)
    pretty_table(io,
        Any[R.θ rad2deg.(R.θ) R.i R.j R.oᵢ R.oⱼ],
        ["Angle (rad)", "Angle (deg)", "Index A", "Index B", "Order A", "Order B"])
    print(io, " case number: $(R.case)\n")
    print(io, " junction:\n   x = $(R.junction[1])\n   y = $(R.junction[2])\n")
end

function BranchingAngleResult(θ::Vector{Float64},
                              i::Vector{Int64},
                              j::Vector{Int64},
                              oᵢ::Vector{Int64},
                              oⱼ::Vector{Int64},
                              junction::NTuple{2,Float64},
                              case::Int64)
    @assert length(θ) == length(i) == length(j) == length(oᵢ) == length(oⱼ)
    BranchingAngleResult(θ, i, j, oᵢ, oⱼ, junction, case, length(θ))
end

function BranchingAngleResult(θ::Float64,
                              i::Int64,
                              j::Int64,
                              oᵢ::Int64,
                              oⱼ::Int64,
                              junction::NTuple{2,Float64},
                              case::Int64)
    BranchingAngleResult([θ], [i], [j], [oᵢ], [oⱼ], junction, case)
end

function BranchingAngleResult(junction::NTuple{2,Float64}, case::Int64)
    BranchingAngleResult([], [], [], [], [], junction, case, 0)
end

Base.isempty(R::BranchingAngleResult) = (R.N == 0)

#-----

function branchinganglecases(geoms, junction, orders, indices)::BranchingAngleResult
    
    #setup
    L = length(geoms)
    idx = sortperm(orders)
    G = geoms[idx]
    O = orders[idx]
    I = indices[idx]
    J = junction

    #gather information this group of stream orders
    omin = minimum(orders)
    omax = maximum(orders)
    nmin = count(o->o==omin, orders)
    nmax = count(o->o==omax, orders)

    #==
    The specific angle calculation(s) depend on the number of streams and their
    orders. The following cases were accumulated to handle the most likely stream
    configurations, but it's likely not perfect or complete. The case numbers 
    used below are leftover from the original Python implementation (not used anymore).
    ==#
    if (L == 2) & (omin == omax == 1)
        #two lonely order 1 streams
        case = 1
        θ = valleyintersectionangle(J, G[1], G[2])
        return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
    elseif (L == 2) & (omin != omax)
        #two lonely streams, one higher order, one lower
        case = 2
        θ = valleyintersectionangle(J, G[1], G[2])
        return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
    elseif (L == 2) & (omin == omax) & (omin != 1)
        #two lonely higher order streams
        case = 3
        θ = valleyintersectionangle(J, G[1], G[2])
        return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
    elseif (L == 3) & (nmin == 2) & (nmax == 1)
        #two lower order flowing into a single higher order
        case = 4
        θ = valleyintersectionangle(J, G[1], G[2])
        return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
    elseif (nmin > 2) & (nmax == 1) & (omin == omax - 1)
        #more than two flowing into a single higher order
        case = 5
        #first calculate the angles of each valley (not intersection angles yet)
        ϕ = [valleyangle(g, J) for g ∈ G]
        #fully rotate angles of low order streams so that their values are greater than the high order angle
        for i ∈ 1:L-1
            @inbounds if ϕ[i] < ϕ[L]
                ϕ[i] += 2π
            end
        end
        #now subtract the higher order stream's angle
        ϕ = ϕ[1:L-1] .- ϕ[L]
        #sort the angles of the lower order streams so that they can be compared in order
        idx = sortperm(ϕ)
        ϕ = ϕ[idx]
        I = I[idx]
        N = length(ϕ)
        #now the intersection angles are in sequence
        θ = ϕ |> diff .|> minorangle .|> abs
        return BranchingAngleResult(θ, I[1:N-1], I[2:N], O[1:N-1], O[2:N], J, case)
    elseif (L == 3) & (nmin == 1) & (nmax == 2)# & (omin == omax -1)
        #single lower order flowing into two higher order (tributary meeting main channel)
        case = 6
        #take the smaller one
        θ₁ = valleyintersectionangle(J, G[1], G[2])
        θ₂ = valleyintersectionangle(J, G[1], G[3])
        if θ₁ < θ₂
            return BranchingAngleResult(θ₁, I[1], I[2], O[1], O[2], J, case)
        else
            return BranchingAngleResult(θ₂, I[1], I[3], O[1], O[3], J, case)
        end
    elseif (L > 2) & (omin == omax) & (omin != 1)
        #more than two higher order, alone?
        case = 7
        return BranchingAngleResult(J, case)
    elseif (L == 4) & (nmin == 2) & (nmax == 2)
        #two tributaries meeting main channel
        case = 8
        #first pair
        θ₁ = [
            valleyintersectionangle(J, G[1], G[3]),
            valleyintersectionangle(J, G[2], G[3])
        ]
        #second pair   
        θ₂ = [
            valleyintersectionangle(J, G[1], G[4]),
            valleyintersectionangle(J, G[2], G[4])
        ]
        #try to take the upstream angles
        if sum(θ₁) < sum(θ₂)
            return BranchingAngleResult(
                θ₁,
                [I[1], I[2]],
                [I[3], I[3]],
                [O[1], O[2]],
                [O[3], O[3]],
                J,
                case
            )
        else
            return BranchingAngleResult(
                θ₂,
                [I[1], I[2]],
                [I[4], I[4]],
                [O[1], O[2]],
                [O[4], O[4]],
                J,
                case
            )           
        end
    elseif (L > 2) & (nmin == 1) & (omax > omin)
        #order 1 joining a variety of higher orders
        case = 9
        #take smallest angle made with next higher order
        idx = findall(o->o==O[2], O)
        θ = Inf
        p = (-1,-1)
        o = (-1, -1)
        for i ∈ idx
            θᵢ = valleyintersectionangle(J, G[1], G[i])
            if θᵢ < θ
                θ = θᵢ
                p = I[1], I[i]
                o = O[1], O[i]
            end
        end
        return BranchingAngleResult(θ, p[1], p[2], o[1], o[2], J, case)
    end

    #empty result
    return BranchingAngleResult(junction, 0)
end

function branchingangles(NA::NetworkAssembly)::Vector{BranchingAngleResult}

    #unpack shapes and junction information
    @unpack geoms, orders, junctions, neighbors = NA
    #number of junctions is same as number of calls to other BranchingAngles() method
    L = length(junctions)
    
    #compute branching angles for each junction asynchronously
    tasks = Vector{Task}(undef, L)
    for i ∈ 1:L
        idx = neighbors[i]
        tasks[i] = @spawn branchinganglecases(geoms[idx], junctions[i], orders[idx], idx)
    end

    #fetch and return
    return BranchingAngleResult[fetch(task) for task ∈ tasks]
end

#--------------------------------------
export dropcases

flatten(x) = collect(Iterators.flatten(x))

function extract(R::Vector{BranchingAngleResult}, f::Symbol)
    flatten(map(r -> getfield(r,f), R))
end

function DataFrames.DataFrame(R::Vector{BranchingAngleResult})
    df = DataFrame(
        "angle" => extract(R, :θ),
        "index A" => extract(R, :i),
        "index B" => extract(R, :j),
        "order A" => extract(R, :oᵢ),
        "order B" => extract(R, :oⱼ),
    )
    df[!,:case] = fill(-1, size(df,1))
    n = 1
    for r ∈ R
        for i ∈ 1:length(r)
            df[n,:case] = r.case
            n += 1
        end
    end
    return df
end

function dropcases(df::DataFrame, cases=Vector{Int64})
    idx = ones(Bool, size(df,1))
    for case ∈ cases
        @. idx &= (df[!,:case] != case)
    end
    return df[idx,:]
end

end
