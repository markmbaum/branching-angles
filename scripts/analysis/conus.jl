using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using CSV
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

df = datadir("exp_pro", "conus_angles.csv") |> CSV.File |> DataFrame
#computes logslope, minorder, maxorder
derivedcols!(df)
#renames precip and temperature columns
renamePT!(df)

##

predictors = [:angle, :logslope, :EVI, :maxorder, :minorder]
targets = [:P, :T, :AI, :SSM, :SUSM]
regions = [:huc2, :huc4, :huc8]

#take a subset of the table
df = df[:,vcat(predictors, targets, regions)]
#standardize all the numeric quantities, leaving the angle unstandardized
transform!(
    df,
    predictors[1:3].=>zscore.=>predictors[1:3],
    targets.=>zscore.=>targets
)

## basic 2D histograms

fig, axs = plt.subplots(
    length(targets),
    length(predictors),
    #sharex=true,
    #sharey=true,
    constrained_layout=true
)
for (i,target) ∈ enumerate(targets)
    for (j,predictor) ∈ enumerate(predictors[1:3])
        Seaborn.histplot(x=df[!,predictor], y=df[!,target], ax=axs[i,j], cmap="magma")
    end
    for (j,predictor) ∈ zip(4:5,predictors[4:5])
        Seaborn.boxplot(x=df[!,predictor], y=df[!,target], ax=axs[i,j])
    end
    axs[i,1].set_ylabel(target)
end
for (j,predictor) ∈ enumerate(predictors)
    axs[end,j].set_xlabel(predictor)
end
foreach(axs) do ax
    ax.set_ylim(-3, 4)
    ax.grid(false)
end

##

hcols = cols[1:end-2]
fig, axs = plt.subplots(2, (length(cols)-2) ÷ 2, sharex=true, sharey=true)
for (ax,col) ∈ zip(axs, hcols)
    Seaborn.histplot(x=df.angle, y=df[!,col], ax=ax, cmap="magma")
    ax.set_title(col)
end
foreach(axs) do ax
    ax.set_ylim(-4, 5)
    ax.grid(false)
end
fig.supxlabel("Branching Angle [radians]")
fig.supylabel("Standardized Environmental Features [-]")
fig.tight_layout()
fig.savefig(plotsdir("histograms"), dpi=500)

## binned average plot & model

m = binstat(df, :angle, cols, mean, 16)
s = binstat(df, :angle, cols, sem, 16)
fig, ax = plt.subplots(1, 1)
for (i,col) in enumerate(cols)
    ax.errorbar(
        m[!,:angle],
        m[!,col],
        s[!,col],
        label=col,
        linewidth=1.5
    )
end
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.82, box.height])
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax.set_xlabel("Bin Average Branching Angle [radians]")
ax.set_ylabel("Bin Average Standardized Response Variables")
fig.savefig(plotsdir("bin_means"), dpi=500)

## correlation matrix

ccols = vcat(:angle,cols)
figure()
ax = Seaborn.heatmap(
    (df[!,ccols] |> Matrix |> cor) .- NaN*I(length(ccols)),
    cmap="RdBu",
    vmin=-1,
    vmax=1
)
ax.set_xticklabels(ccols, rotation=45)
ax.set_yticklabels(ccols, rotation=0)
ax.set_title("Pearson Correlation Matrix")
ax.figure.tight_layout()
ax.figure.savefig(plotsdir("corr_matrix"), dpi=500)

## standardize the angle for regression

transform!(df, :angle=>zscore=>:angle)

## linear models

ecols = [:angle, :logslope, :EVI, :maxorder, :minorder]
pcols = [:P, :AI, :SSM, :SUSM]
lms = map(pcols) do col
    lm(term(col) ~ term(0) + (map(term, ecols) |> sum), df)
end

fig, ax = plt.subplots(1,1)
x = collect(1:length(ecols))
w = 0.15
for (i,model) ∈ enumerate(lms)
    ax.bar(x .+ w*(i-1), coef(model), width=w, label=pcols[i])
end
ax.legend()
ax.set_xticks(x .+ 0.15*4/2)
ax.set_xticklabels(ecols)
ax.grid(axis="x", false)
ax.set_ylabel("Linear Regression Coefficient")
fig.savefig(plotsdir("linear_models"), dpi=500)

##

ecols = [:angle, :logslope, :maxorder, :minorder]
tcols = term.(ecols)
pcols = [:P, :EVI, :SSM, :SUSM, :AI]
lms = map(pcols) do col
    lm(term(col) ~ sum(tcols) + sum(tcols .& tcols), df)
end

fig, axs = plt.subplots(
    length(pcols),
    length(regions),
    constrained_layout=true,
    sharex=true
)
for (j,region) ∈ enumerate(regions)
    g = combine(
        groupby(df, region),
        pcols.=>mean.=>pcols,
        ecols.=>mean.=>ecols,
        nrow
    )
    filter!(r -> r.nrow > 25, g)
    println(median(g.nrow))
    for (i,col) ∈ enumerate(pcols)
        p = predict(lms[i], g[!,ecols])
        r = g[!,col] .- p
        println("  $(std(r))")
        Seaborn.histplot(r, ax=axs[i,j])
        axs[i,j].set_ylabel("")
    end
end
for (i,col) ∈ enumerate(pcols)
    axs[i,1].set_ylabel(col)
end
for (i,region) ∈ enumerate(regions)
    axs[1,i].set_title(region)
end

#fig, axs = plt.subplots(1, length(regions))
#for i ∈ 1:length(regions)
#    g = combine(groupby(df, region), nrows)
#    Seaborn.histplot(g.nrows, ax=axs[i])
#    
#end

## lasso paths

ecols = [:angle, :logslope, :maxorder, :minorder]
pcols = [:P, :EVI, :SSM, :SUSM, :AI]
lps = map(pcols) do col
    fit(
        LassoPath,
        term(col) ~ map(term, ecols) |> sum,
        df,
        intercept=false
    )
end

fig, axs = plt.subplots(1, length(pcols), figsize=(10,4), constrained_layout=true)
v = map(x -> x |> coef .|> abs |> maximum, lps) |> maximum
for i ∈ 1:length(pcols)
    c = axs[i].pcolormesh(
        lps[i].model.λ,
        1:length(ecols),
        coef(lps[i]),
        vmin=-v,
        vmax=v,
        cmap="RdBu"
    )
    axs[i].set_title(pcols[i])
    axs[i].set_yticks(1:length(ecols))
end
axs[1].set_yticklabels(ecols)
foreach(ax->ax.set_yticklabels([]), axs[2:end])
cb = plt.colorbar(c, ax=axs)
cb.set_label("Lasso Coefficient", rotation=270, va="bottom")
fig.supxlabel("Regularization Parameter (λ)")
