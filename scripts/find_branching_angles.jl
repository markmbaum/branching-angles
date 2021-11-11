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

#=======================
MARTIAN BRANCHING ANGLES
=======================#

## information about which streams are at which junctions

fn = datadir("exp_raw", "mars_networks", "Hynek_Valleys_geology_draped_mercator.shp")
marsassembly = assemblenetworks(fn);
fn = datadir("exp_pro", "mars_network_assembly")
serialize(fn, marsassembly)

## all the branching angles as a vector of one BranchingAngleResult struct for each junction

marsangles = branchingangles(marsassembly, :VNOrder)
fn = datadir("dir_pro", "mars_branching_angles")
serialize(fn, marsangles)

#========================================
CONTIGUOUS UNITED STATES BRANCHING ANGLES
========================================#

## information about which streams are at which junctions
fn = datadir("exp_raw", "conus-networks", "NHDFlowline_Network_NoMZ_w_basins.shp")
conusassembly = assemblenetworks(fn)
fn = datadir("dir_pro", "conus_network_assembly")
serialize(fn, conusassembly)

## all the branching angles as a vector of one BranchingAngleResult struct for each junction
conusangles = branchingangles(conusassembly, :StreamOrde)
fn = datadir("dir_pro", "conus_branching_angles")
serialize(fn, conusangles)