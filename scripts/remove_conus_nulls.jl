using DrWatson
@quickactivate "Martian Branching Angles"
using DataFrames
using CSV

##

fname = datadir("exp_pro", "conus_angles.csv")

## load the nearly finalized table

df = CSV.read(fname, DataFrame)

## build up mask representing rows with nulls

b = (
    @. (df[!,:ppt_annual] == -100.0) |
    #NDVI and EVI should have NaN in the same places
    isnan(df[!,:NDVI]) |
    #SMAP columns
    isnan(df[!,:SSM]) | isnan(df[!,:SUSM]) | isnan(df[!,:SMP]) |
    #the state boundaries missed a few
    ismissing(df[!,:state]) |
    #there is a small number of null slopes
    isnan(df[!,"slope A"]) | isnan(df[!,"slope B"]) |
    #a few nulls in the aridity index
    ismissing(df[!,:AI])
)

println("$(sum(b)) nulls out of $(length(b)) rows")

## remove the rows and save with a new file name

CSV.write(replace(fname, ".csv"=>"_nulls_removed.csv"), df[Not(b),:])