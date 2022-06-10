using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Arrow
using DataFrames
using GLM
using Lasso
using LinearAlgebra: I
using PyPlot
using Seaborn
using Statistics
using StatsBase

pygui(true)

## load and prepare CONUS data

df = datadir("exp_pro", "conus_angles.feather") |> Arrow.Table |> DataFrame
df = mapcols(Vector, df)
#computes logslope, minorder, maxorder
derivedcols!(df)
#renames precip and temperature columns
renamePT!(df)

##

fig, ax = plt.subplots(1, 1, figsize=(7,3), constrained_layout=true)
Seaborn.histplot(x=df.angle, ax=ax, color="gray", stat="count")
ax.set_xlim(0, π)
ax.set_xlabel("Branching Angle [radians]")
fig.savefig(plotsdir("angle_histogram"), dpi=500)

##

pred = [:P, :T, :AI, :EVI, :SSM, :logslope, :minorder, :maxorder]
cont = pred[1:6]
regions = [:huc2, :huc4, :huc8]

select!(df, vcat(:angle, pred, regions))
transform!(
    df,
    cont.=>zscore.=>cont,
    :angle=>zscore=>:angle
)
filter!(r -> r.minorder < 7, df) #there are too few instances with minorder >= 7

## 2D histograms

fig, axs = plt.subplots(
    1, 6,
    figsize=(9,4),
    sharex=true,
    sharey=true, 
    constrained_layout=true
)
for (ax,c) ∈ zip(axs, cont)
    Seaborn.histplot(x=df.angle, y=df[!,c], cmap="magma", ax=ax)
    ax.set_title(c)
end
axs[1].set_ylim(-4, 4)
fig.supxlabel("Branching Angle")
fig.supylabel("Standardized Predictors")
fig.savefig(plotsdir("histograms"), dpi=500)

## violins for stream orders

fig, axs = plt.subplots(
    1, 2, 
    figsize=(7,5), 
    sharex=true, 
    sharey=true, 
    constrained_layout=true
)
sl = df[df.AI .<= -1,:]
for (ax,c) ∈ zip(axs, [:minorder, :maxorder])
    Seaborn.violinplot(
        x=sl.angle,
        y=sl[!,c],
        ax=ax,
        orient="h",
        #scale="count",
        palette="crest"
    )
    ax.set_title(c)
end
axs[1].set_title("Minimum Stream Order")
axs[2].set_title("Maximum Stream Order")
fig.supxlabel("Branching Angle")
fig.savefig(plotsdir("violins"), dpi=500)

## binned average plot & model

m = binstat(mapcols(zscore, df), :angle, pred, mean, 16)
s = binstat(mapcols(zscore, df), :angle, pred, sem, 16)
fig, ax = plt.subplots(1, 1, figsize=(7,4), constrained_layout=true)
for (i,c) in enumerate(pred)
    ax.errorbar(
        m[!,:angle],
        m[!,c],
        s[!,c],
        label=c,
        linewidth=1.5
    )
end
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.82, box.height])
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax.set_xlabel("Bin Average Standardized Branching Angle [radians]")
ax.set_ylabel("Bin Average Standardized Response Variables")
fig.savefig(plotsdir("bin_means"), dpi=500)

## correlation matrix

c = vcat(:angle,pred)
figure()
M = (df[!,c] |> Matrix |> cor)
for i in 1:size(M,1), j in i:size(M,2)
    M[i,j] = NaN
end
M = M[2:end,1:end-1]
ax = Seaborn.heatmap(
    M,
    cmap="RdBu",
    vmin=-1,
    vmax=1
)
ax.set_xticklabels(c[1:end-1], rotation=45)
ax.set_yticklabels(c[2:end], rotation=0)
ax.set_title("Pearson Correlation Matrix")
ax.figure.tight_layout()
ax.figure.savefig(plotsdir("corr_matrix"), dpi=500)

## bin by watershed and take the means

ge = combine(
    groupby(
        df,
        :huc8
    ),
    pred .=> mean .=> pred,
    :angle => mean => :angle,
    nrow => :n
)
filter!(r -> r.n > 10, ge)
transform!(ge, pred .=> zscore .=> pred)

## huc8 histograms

fig, axs = plt.subplots(
    2, 4,
    figsize=(9,4),
    sharex=true,
    sharey=true,
    constrained_layout=true
)
for (ax,c) ∈ zip(axs, pred)
    Seaborn.histplot(x=ge.angle, y=ge[!,c], ax=ax, cmap="magma")
    Seaborn.regplot(
        x=ge.angle,
        y=ge[!,c],
        ax=ax,
        robust=true,
        scatter=false,
        ci=nothing,
        color="deepskyblue"
    )
    ax.set_title(c)
end
axs[1].set_ylim(-4, 4)
fig.supxlabel("Standardized Branching Angle")
fig.supylabel("Standardized Predictors")
fig.savefig(plotsdir("histograms_huc8"), dpi=500)

## lasso paths

ldf = ge[:,vcat(pred, :angle)]
for c in string.(pred)
    ldf[!,c*"²"] = ldf[!,c].^2
end
transform!(ldf, names(ldf) .=> zscore .=> names(ldf))
cols = filter(x -> x != "angle", names(ldf))
lp = StatsModels.fit(
    LassoPath,
    ldf[:,cols] |> Matrix,
    ldf.angle,
    intercept=false
)
c = Lasso.coef(lp) .|> abs

fig, ax = plt.subplots(1, 1, constrained_layout=true)
r = ax.pcolormesh(
    lp.λ,
    1:size(c,1),
    c,
    vmin=0,
    cmap="magma"
)
ax.set_title("Lasso Path Coefficient Magnitudes")
plt.colorbar(r, ax=ax)
ax.set_xscale("log")
ax.set_xlabel("Regularization Parameter (λ)")
ax.set_yticks(1:size(c,1))
ax.set_yticklabels(cols)
#fig.savefig(plotsdir("lasso_path"), dpi=500)

## lasso paths

ldf = df
uminorder = ldf.minorder |> unique |> sort
umaxorder = ldf.maxorder |> unique |> sort
lp = fit(
    LassoPath,
    hcat(
        ldf[!,cont] |> Matrix,
        (uminorder .== permutedims(ldf.minorder))',
        (umaxorder .== permutedims(ldf.maxorder))'
    ),
    ldf.angle,
    intercept=false
)
c = coef(lp) .|> abs

fig, ax = plt.subplots(1, 1, constrained_layout=true)
r = ax.pcolormesh(
    lp.λ,
    1:size(c,1),
    c,
    vmin=0,
    cmap="magma"
)
ax.set_title("Lasso Path Coefficient Magnitudes")
plt.colorbar(r, ax=ax)
ax.set_xscale("log")
ax.set_xlabel("Regularization Parameter (λ)")
ax.set_yticks(1:size(c,1))
ax.set_yticklabels(vcat(
    cont,
    ["minorder "*string(o) for o in uminorder],
    ["maxorder "*string(o) for o in umaxorder]
))
fig.savefig(plotsdir("lasso_path"), dpi=500)