using DrWatson
@quickactivate "Martian Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Serialization

##

#===
The network assembly should already be done by the assemble_networks.jl
script, which is memory and time intensive, so done on the cluster.
The branching angles are relatively quick to compute once the network
information has been assembled, so it's done in this script, which can
be run locally on the laptop. Branching angles are computed in parallel,
asynchronously, which helps speed things up but requires Julia to be
launched with multiple threads to take advantage.
===#

# fnassembly - path to serialized NetworkAssembly struct
# fnanbles - output path for serialized vector of BranchingAngleResult structs
function wrangleangles(fnassembly, fnangles)::Nothing
    #one fell swoop
    serialize(fnangles, branchingangles(deserialize(fnassembly)))
    println("file created: $fnangles")
    return nothing
end

## Mars angles

@time wrangleangles(
    datadir("exp_pro", "mars_assembly"),
    datadir("exp_pro", "mars_angles")
)

## CONUS angles

@time wrangleangles(
    datadir("exp_pro", "conus_assembly"),
    datadir("exp_pro", "conus_angles")
)