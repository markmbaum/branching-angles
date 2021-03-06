module BranchingAngles

using Base.Threads: @threads, Task, @spawn, fetch
using Base.Iterators: partition
using CSV
using DataFrames: DataFrame, transform!, rename!, groupby, combine
import DataFrames
using DrWatson
using Graphs
using MultiAssign
using Optim
using PrettyTables
using Shapefile: Table
using StatsBase
using UnPack

#------------------------------------------------------------------------------
export endintersection, endintersecting
export assemblenetworks
export NetworkAssembly

struct NetworkAssembly{T}
    geoms::Vector{T} #shapefile geometries
    orders::Vector{Int64} #stream order of each geometry
    networks::Vector{Vector{Int64}} #groups of stream/valley indices for each full network
    junctions::Vector{NTuple{2,Float64}} #all junction points
    neighbors::Vector{Vector{Int64}} #stream/valley indices meeting at each junction
end

function Base.show(io::IO, A::NetworkAssembly{T}) where {T}
    println(io, "NetworkAssembly:")
    println(io, "  geometries of type $T")
    println(io, "  $(length(A.junctions)) junctions")
    println(io, "  $(length(A.networks)) networks")
    print(io, "  stream orders $(sort(unique(A.orders)))")
end

function endintersecting(a₁, aₙ, b₁, bₙ)::Bool
    if (a₁ == b₁) || (a₁ == bₙ) || (aₙ == b₁) || (aₙ == bₙ)
        return true
    end
    return false
end

function endintersection(a₁, aₙ, b₁, bₙ)
    if a₁ == b₁
        return a₁
    elseif a₁ == bₙ
        return a₁
    elseif aₙ == b₁
        return aₙ
    elseif aₙ == bₙ
        return aₙ
    end
    error("no intersecting end point found")
end

function endintersecting(a::AbstractVector{Float64},
                         b::AbstractVector{Float64})::Bool
    @inbounds endintersecting(
        (a[1], a[2]),
        (a[3], a[4]),
        (b[1], b[2]),
        (b[3], b[4])
    )
end

function endintersection(a::AbstractVector{Float64},
                         b::AbstractVector{Float64})::NTuple{2,Float64}
    @inbounds begin
        a₁ = a[1], a[2]
        aₙ = a[3], a[4]
        b₁ = b[1], b[2]
        bₙ = b[3], b[4]
    end
    endintersection(a₁, aₙ, b₁, bₙ)
end

function endpointmatrix(geoms::Vector{T})::Matrix{Float64} where {T}
    #total number of geometries
    N = length(geoms)
    #put all the end points into one array/matrix
    E = zeros(Float64, 4, N)
    @inbounds @threads for i ∈ 1:N
        p = geoms[i].points
        L = length(p)
        E[1,i] = p[1].x
        E[2,i] = p[1].y
        E[3,i] = p[L].x
        E[4,i] = p[L].y
    end
    return E
end

function adjacencylists(E::Matrix{Float64})::Vector{Set{Int64}}
    #number of geometries
    @assert size(E,1) == 4
    N = size(E,2)
    #initalize adjacency lists
    I = [Set{Int64}() for _ ∈ 1:N]
    #generate adjacency lists in parallel (not sure if graphs are thread-safe)
    @inbounds @threads for i ∈ 1:N
        #see if geometries intersect
        gᵢ = @view E[:,i]
        for j ∈ 1:N
            gⱼ = @view E[:,j]
            if endintersecting(gᵢ, gⱼ)
                push!(I[i], j)
            end
        end
    end
    return I
end

function adjacencygraph(I::Vector{Set{Int64}})::SimpleGraph{Int64}
    #create the graph
    N = length(I)
    G = Graph(N)
    @inbounds for i ∈ 1:N
        for j ∈ I[i]
            add_edge!(G, i, j)
        end
    end
    return G
end

function assemblenetworks(E::Matrix{Float64})
    @assert size(E,1) == 4
    #create adjacency graph and dissolve it all at once
    E |> adjacencylists |> adjacencygraph |> connected_components
end

function findjunctions(network::Vector{Int64}, F::AbstractMatrix{Float64})
    @assert size(F,1) == 4
    Lₙ = length(network)
    x = Float64[] #first coordinate of junction points
    y = Float64[] #second coordinate of junction points
    s = Set{NTuple{2,Float64}}() #for preventing duplicates
    neighbors = Vector{Int64}[]
    @inbounds if Lₙ > 1
        #find all junction points in the group
        for j ∈ 1:Lₙ-1
            gⱼ = @view F[:,j]
            for k ∈ j+1:Lₙ
                gₖ = @view F[:,k]
                #see if there is an intersection point
                if endintersecting(gⱼ, gₖ)
                    #extract the intersection point
                    p = endintersection(gⱼ, gₖ)
                    #do not add a duplicate
                    if !(p ∈ s)
                        push!(x, p[1])
                        push!(y, p[2])
                        push!(s, p)
                    end
                end
            end
        end
        #now find which lines meet at the newly identified junctions
        for j ∈ 1:length(x)
            #junction j coordinates
            pⱼ = x[j], y[j]
            #start a new neighborhood around the junction
            push!(neighbors, Int64[])
            #check which geometries are in the hood
            for i ∈ 1:length(network)
                a = F[1,i], F[2,i]
                b = F[3,i], F[4,i]
                if (pⱼ == a) | (pⱼ == b)
                    push!(neighbors[end], network[i])
                end
            end
        end
    end
    return x, y, neighbors
end

function assemblenetworks(geoms::Vector{T}) where {T}
    #end points only
    E = endpointmatrix(geoms)
    #get network groups
    networks = assemblenetworks(E)
    N = length(networks)
    println(stdout, "$N networks identified")
    #find junction points and valleys meeting at each junction
    tasks = Vector{Task}(undef,N)
    for i ∈ 1:N
        tasks[i] = @spawn findjunctions(networks[i], E[:,networks[i]])
    end
    results = map(fetch, tasks)
    #structure the results
    x = flatten(map(r->r[1], results))
    y = flatten(map(r->r[2], results))
    junctions = collect(zip(x,y))
    neighbors = flatten(map(r->r[3], results))
    println(stdout, "$(length(junctions)) total junctions identified")
    flush(stdout)
    return networks, junctions, neighbors
end

function assemblenetworks(geoms::Vector{T}, orders::Vector{Int64})::NetworkAssembly where {T}
    @assert length(geoms) == length(orders)
    networks, junctions, neighbors = assemblenetworks(geoms)
    #final construction
    NetworkAssembly(
        geoms,
        orders,
        networks,
        junctions,
        neighbors
    )
end

function assemblenetworks(fn::String, ordercol::Symbol)::NetworkAssembly
    println(stdout, "assembling networks from file: $fn")
    table =  fn |> Table |> DataFrame
    geoms = table[:,:geometry]
    orders = table[:,ordercol]
    println(stdout, "$(length(geoms)) geometries present")
    println(stdout, "stream orders present: $(sort(unique(orders)))")
    flush(stdout)
    assemblenetworks(geoms, orders)
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
(vc::ValleyCoordinates)(m, b) = orthogonalloss(m, b, vc.x, vc.y, vc.N)
(vc::ValleyCoordinates)(X) = vc(X[1], X[2])

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

getx(geom) = [p.x for p ∈ geom.points]
gety(geom) = [p.y for p ∈ geom.points]

valleyslope(geom) = valleyslope(getx(geom), gety(geom), geom.measures)

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

function valleyangle(geom, J::NTuple{2,Float64})
    valleyangle(getx(geom), gety(geom), J[1], J[2])
end

function minorangle(θ)
    while (θ > π) | (θ < -π)
        θ = (θ > π) ? θ - 2π : θ + 2π
    end
    return θ
end

function valleyintersectionangle(x₁::AbstractVector, y₁::AbstractVector, J₁::NTuple{2,Float64},
                                 x₂::AbstractVector, y₂::AbstractVector, J₂::NTuple{2,Float64})
    #direction angle of each valley
    θ₁ = valleyangle(x₁, y₁, J₁[1], J₁[2])
    θ₂ = valleyangle(x₂, y₂, J₂[1], J₂[2])
    #difference between the angles
    ϕ = θ₁ - θ₂
    #correct for results θ > π or θ < -π
    return abs(minorangle(ϕ))
end

function valleyintersectionangle(J, geom₁, geom₂)
    valleyintersectionangle(
        getx(geom₁), gety(geom₁), J,
        getx(geom₂), gety(geom₂), J
    )
end

#-----

struct BranchingAngleResult
    θ::Vector{Float64} #branching angles
    i::Vector{Int64} #indices of first shape in pairs
    j::Vector{Int64} #indices of second shape in pairs
    oᵢ::Vector{Int64} #stream orders of shapes i
    oⱼ::Vector{Int64} #stream orders of shapes j
    junction::NTuple{2,Float64} #junction location
    case::Int64 #branching case/type
    N::Int64 #vector lengths (number of angles)
end

Base.length(R::BranchingAngleResult) = R.N

function Base.show(io::IO, R::BranchingAngleResult)
    pretty_table(io,
        Any[R.θ rad2deg.(R.θ) R.i R.j R.oᵢ R.oⱼ],
        ["Angle (rad)", "Angle (deg)", "Index A", "Index B", "Order A", "Order B"])
    print(io, " case number: $(R.case)\n")
    print(io, " junction:\n   x = $(R.junction.x)\n   y = $(R.junction.y)\n")
    print(io, " $(R.njunc) geometries at junction")
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
    @assert length(geoms) == length(indices) == length(orders)
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
    if L == 2
        if omin == omax == 1
            #two lonely order 1 streams
            case = 1
            θ = valleyintersectionangle(J, G[1], G[2])
            return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
        elseif omin != omax
            #two lonely streams, one higher order, one lower
            case = 2
            θ = valleyintersectionangle(J, G[1], G[2])
            return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
        elseif (omin == omax) & (omin != 1)
            #two lonely higher order streams
            case = 3
            θ = valleyintersectionangle(J, G[1], G[2])
            return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
        end
    elseif L == 3
        if nmin == 2
            #two lower order flowing into a single higher order
            case = 4
            θ = valleyintersectionangle(J, G[1], G[2])
            return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
        elseif (nmin == 1) & (nmax == 2)# & (omin == omax -1)
            #single lower order flowing into two higher order (tributary meeting main channel)
            case = 5
            #take the smaller one
            θ₁ = valleyintersectionangle(J, G[1], G[2])
            θ₂ = valleyintersectionangle(J, G[1], G[3])
            if θ₁ < θ₂
                return BranchingAngleResult(θ₁, I[1], I[2], O[1], O[2], J, case)
            else
                return BranchingAngleResult(θ₂, I[1], I[3], O[1], O[3], J, case)
            end
        elseif (nmin == 1) & (nmax == 1)
            #three of different orders
            case = 6
            #take the first pair
            θ = valleyintersectionangle(J, G[1], G[2])
            return BranchingAngleResult(θ, I[1], I[2], O[1], O[2], J, case)
        end
    elseif (nmin > 2) & (nmax == 1) & (omin == omax - 1)
        #more than two flowing into a single higher order
        case = 7
        #first calculate the angles of each valley (not intersection angles yet)
        ϕ = [valleyangle(g, J) for g ∈ G]
        #fully rotate angles of low order streams so that their values are greater than the high order angle
        for i ∈ 1:L-1
            if ϕ[i] < ϕ[L]
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
    elseif (L > 2) & (omin == omax) & (omin != 1)
        #more than two higher order, alone?
        case = 10
        return BranchingAngleResult(J, case)
    end

    return BranchingAngleResult(J, 0)
end

function branchingangles(NA::NetworkAssembly)::Vector{BranchingAngleResult}

    #unpack shapes and junction information
    @unpack geoms, orders, junctions, neighbors = NA
    #check each junction for branching angles
    Lⱼ = length(junctions)
    
    #compute branching angles for each junction asynchronously
    res = Vector{BranchingAngleResult}(undef,Lⱼ)
    @threads for i ∈ 1:Lⱼ
        J = junctions[i]
        n = neighbors[i]
        res[i] = branchinganglecases(geoms[n], J, orders[n], n)
    end
    return res
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
    L = size(df,1)
    df[!,:case] = fill(-1, L)
    @multiassign df[!,:x], df[!,:y] = fill(NaN, L)
    i = 1
    for r ∈ R
        for _ ∈ 1:r.N
            df[i,:case] = r.case
            df[i,:x] = r.junction[1]
            df[i,:y] = r.junction[2]
            i += 1
        end
    end
    return df
end

function dropcases(df::DataFrame, cases::Int...)
    idx = ones(Bool, size(df,1))
    for case ∈ cases
        @. idx &= (df[!,:case] != case)
    end
    return df[idx,:]
end

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# functions for preparing final data

export logslope!, maporder!, derivedcols!, standardize!, renamePT!, binstat

function logslope!(df::DataFrame)::Nothing
    a = df[:,"slope A"]
    b = df[:,"slope B"]
    a[a .< 1e-5] .= 1e-5
    b[b .< 1e-5] .= 1e-5
    df[!,:logslope] = @. log10(a)/2 + log10(b)/2
    nothing
end

function maporder!(df::DataFrame, f::F)::Nothing where {F<:Function}
    df[!,string(f)*"order"] = map(x->f(x...), zip(df[!,"order A"], df[!,"order B"]))
    nothing
end

function derivedcols!(df::DataFrame)::Nothing
    logslope!(df)
    maporder!(df, max)
    maporder!(df, min)
    nothing
end

function renamePT!(df::DataFrame)::Nothing
    for col ∈ names(df)
        if occursin("ppt", col)
            rename!(df, col => replace(col, "ppt"=>"P", "_annual"=>"", "_"=>""))
        elseif occursin("tmean", col)
            rename!(df, col => replace(col, "tmean"=>"T", "_annual"=>"", "_"=>""))
        end
    end
    nothing
end

function binstat(df::DataFrame,
                 bincol::Union{Symbol,String},
                 statcols::AbstractArray,
                 stat::Function,
                 nbins::Int=10)
    df = copy(df)
    h = fit(Histogram, df[!,bincol], nbins=nbins)
    df[!,:binindex] = map(x->StatsBase.binindex(h,x), df[!,bincol])
    combine(
        groupby(
            df,
            :binindex
        ),
        bincol => stat => bincol,
        statcols .=> stat .=> statcols
    )
end

end