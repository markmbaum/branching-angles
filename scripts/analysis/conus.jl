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
#SMP is nearly perfectly correlated with SSM, only need one of them
select!(df, Not(:SMP))

## a single huge lasso path

dg = select(df, r"^(angle|P|T|AI|EVI|SSM|logslope|minorder|maxorder)")
mapcols!(zscore, dg)
predictors = names(dg)[2:end]

lp = fit(
    LassoPath,
    term(:angle) ~ sum(term.(predictors)),
    dg,
    intercept=false,
    nλ=100
)

fig = figure(figsize=(5,7))
ax = Seaborn.heatmap(abs.(coef(lp)), cmap="magma")

ax.set_title("Lasso Coefficient Magnitude")
ax.set_xlabel("Regularization Parameter (λ)")
λ = round.(lp.model.λ, sigdigits=3)
ax.set_xticks([0,length(λ)-1])
ax.set_xticklabels([λ[1], λ[end]], rotation=45)

ax.set_yticks(collect(0:length(predictors)-1) .+ 0.5)
ax.set_yticklabels(predictors, rotation=0)
ax.set_ylabel("feature/predictor")

fig.tight_layout()
fig.savefig(plotsdir("lasso_path_big"), dpi=500)

##

#state and watershed labels
regions = [
    :huc4,
    :huc6,
    :huc8
]
#columns of primary interest
cols = [
    :T,
    :P,
    :EVI,
    :SSM,
    :SUSM,
    :AI,
    :logslope,
    :maxorder,
    :minorder
]
#take a subset of the table
df = df[:,vcat(:angle, cols, regions)]
#standardize all the numeric quantities
transform!(df, cols[1:end-2].=>zscore.=>cols[1:end-2])

## order distributions

sf = stack(
    filter(r -> r.minorder < 7, df), #only 2 minorder == 7
    [:minorder, :maxorder],
    :angle
)
ax = Seaborn.lineplot(
    x=sf.value,
    y=sf.angle,
    hue=sf.variable
)
ax.set_xlabel("Stream Order Metric at Junction")
ax.set_ylabel("Mean Branching Angle")
ax.figure.savefig(plotsdir("order_metrics"), dpi=500)

## standardize the order metrics

transform!(df, [:minorder,:maxorder].=>zscore.=>[:minorder,:maxorder])

## basic histograms

hcols = cols[1:end-2]
fig, axs = plt.subplots(2, Int(ceil(length(hcols)/2)), sharex=true, sharey=true)
for (ax,col) ∈ zip(axs, hcols)
    Seaborn.histplot(x=df.angle, y=df[!,col], ax=ax, cmap="magma")
    ax.set_title(col)
end
foreach(ax->ax.set_ylim(-4, 5), axs)
fig.supxlabel("Branching Angle [radians]")
fig.supylabel("Standardized Environmental Features [-]")
fig.tight_layout()
axs[end].set_visible(false)
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
ax = Seaborn.heatmap(
    (df[!,vcat(ccols)] |> Matrix |> cor) .- NaN*I(length(ccols)),
    cmap="RdBu_r",
    vmin=-1,
    vmax=1
)
ax.set_xticklabels(ccols, rotation=45)
ax.set_yticklabels(ccols, rotation=0)
ax.set_title("Pearson Correlation Matrix")
ax.figure.tight_layout()
ax.figure.savefig(plotsdir("corr_matrix"), dpi=500)

## multiple linear regression and regional robustness tests

rcols = cols[:]
tcols = term.(rcols)
model = fit(
    LinearModel,
    term(:angle) ~ sum(tcols),
    df
)
Base.show(model)

fig, axs = plt.subplots(1, 2, figsize=(7,4), sharex=true)
for (region, ax) ∈ zip([:huc4, :huc8], vec(axs))
    g = combine(
        groupby(
            df,
            region
        ),
        rcols.=>mean.=>rcols,
        :angle=>mean=>:angle,
        nrow
    )
    println("\nmedian angle: $(median(g.angle))")
    filter!(r -> r.nrow > 25, g)
    p = predict(model, hcat(ones(size(g,1)), g[!,rcols]))
    resid = p .- g.angle
    Seaborn.histplot(x=resid, kde=true, ax=ax)
    ax.set_title(region)
    ax.set_xlabel("")
    ax.set_ylabel("")
end
fig.supxlabel("Branching Angle Residual [radians]")
fig.supylabel("Count")
fig.tight_layout()
fig.savefig(plotsdir("region_residuals"), dpi=500)

## vegetation robustness test

σs = LinRange(0, 1, 20)
ms4 = similar(σs)
ms8 = similar(σs)
for (i,σ) ∈ enumerate(σs)
    model = fit(
        LinearModel,
        term(:angle) ~ sum(term.(cols)),
        filter(r -> r.EVI > σ, df)
    )
    g4 = combine(
        groupby(
            df,
            :huc4
        ),
        cols.=>mean.=>cols,
        :angle=>mean=>:angle,
        nrow
    )
    sl = filter(r -> r.EVI < σ, g4)
    p = predict(model, hcat(ones(size(sl,1)), sl[!,cols]))
    resid = p .- sl.angle
    ms4[i] = mean(resid)
    g8 = combine(
        groupby(
            df,
           :huc8
        ),
        cols.=>mean.=>cols,
        :angle=>mean=>:angle,
        nrow
    )
    
    sl = filter(r -> r.EVI < σ, g8)
    p = predict(model, hcat(ones(size(sl,1)), sl[!,cols]))
    resid = p .- sl.angle
    ms8[i] = mean(resid)
end
plt.figure()
plt.errorbar(σs, ms4, σ, ms8)

##

lcols = cols[:]
tcols = term.(lcols)
lpl = fit(
    LassoPath,
    term(:angle) ~ sum(tcols),
    df[!,vcat(:angle,lcols)],
    intercept=false
)
lpq = fit(
    LassoPath,
    term(:angle) ~ sum(tcols) + sum(tcols .& tcols),
    df[!,vcat(:angle,lcols)],
    intercept=false
)

figure()
axl = Seaborn.heatmap(abs.(coef(lpl)), cmap="magma")
figure()
axq = Seaborn.heatmap(abs.(coef(lpq)), cmap="magma")

for (ax,lp) ∈ zip((axl,axq), (lpl,lpq))
    ax.set_title("Abs. Lasso Coefficients")
    ax.set_xlabel("Regularization Parameter (λ)")
    λ = round.(lp.model.λ, sigdigits=3)
    ax.set_xticks([0,length(λ)-1])
    ax.set_xticklabels([λ[1], λ[end]], rotation=45)
end

axl.set_yticklabels(lcols, rotation=0)
axq.set_yticklabels(vcat(lcols, map(c->string(c)*"²", lcols)), rotation=0)
foreach(ax->ax.figure.tight_layout(), (axl,axq))
axl.figure.savefig(plotsdir("lasso_path_linear"), dpi=500)
axq.figure.savefig(plotsdir("lasso_path_quadratic"), dpi=500)