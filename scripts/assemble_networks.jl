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
directory as serialized files. The same applies to all shapefiles.
===#

#Martian networks
serialize(
	datadir(
		"exp_pro",
		"mars_assembly_serialized"
	),
	assemblenetworks(
		datadir(
			"exp_raw",
			"mars-networks",
			"Hynek_Valleys_geology_draped_mercator.shp"
		),
		:VNOrder
	)
)

#Contiguous United States (CONUS) networks
serialize(
	datadir(
		"exp_pro",
		"conus_assembly_serialized"
	),
	assemblenetworks(
		datadir(
			"exp_raw",
			"conus-networks",
			"NHDFlowline_Network_NoMZ_w_basins.shp"
		),
		:StreamOrde #[sic]
	)
)
