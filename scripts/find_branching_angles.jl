using DrWatson
@quickactivate "Martian Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Serialization

#===
To calculate branching angles, first we need to gather some information 
about the valley networks. There are two critical groups of information
for the branching angle calculations that come next.
  1. All the junction points where valleys intersect
  2. The indices of the valleys intersecting at those junction points
The following code assembles that information and writes it into the "processed data"
directory as text files. The same applies to all shapefiles.
===#

function findbranchingangles(fnshp::String,
						     order::Symbol,
							 fnassembly::String,
							 fnangles::String
							 )::Nothing
	assembly = assemblenetworks(fnshp, order)
	serialize(fnassembly, assembly)
	angles = branchingangles(assembly)
	serialize(fnangles, angles)
	return nothing
end

##

#=======================
MARTIAN BRANCHING ANGLES
=======================#

findbranchingangles(
	datadir(
		"exp_raw",
		"mars-networks",
		"Hynek_Valleys_geology_draped_mercator.shp"
	),
	:VNOrder,
	datadir(
		"exp_pro",
		"mars_assembly_serialized"
	),
	datadir(
		"exp_pro",
		"mars_angles_serialized"
	)
)

##

#========================================
CONTIGUOUS UNITED STATES BRANCHING ANGLES
========================================#

findbranchingangles(
	datadir(
		"exp_raw",
		"conus-networks",
		"NHDFlowline_Network_NoMZ_w_basins.shp"
	),
	:StreamOrde, #[sic]
	datadir(
		"exp_pro",
		"conus_assembly_serialized"
	),
	datadir(
		"exp_pro",
		"conus_angles_serialized"
	)
)

