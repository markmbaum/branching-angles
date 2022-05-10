using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using CSV
using DataFrames
using Flux
using Flux: mse, params, update!, onehotbatch
using Base.Iterators: partition
using PyPlot
using Random: shuffle
using StatsBase
using Statistics

pygui(true)

## FUNCTIONS

function shuffledata(X, Y)
    @assert size(X,2) == size(Y,2)
    idx = shuffle(1:size(X,2))
    X[:,idx], Y[:,idx]
end

function splittenths(X, Y, tenth=9)
    M = size(X,2)
    idx = tenth*Int(M ÷ 10)
    X[:,1:idx], Y[:,1:idx], X[:,idx+1:end], Y[:,idx+1:end]
end

function train(X, Y, model, loss, opt, epochs, nbatch)
    #data shape
    M = size(X,2)
    @assert size(Y,2) == M
    #split for train/validation
    Xtr, Ytr, Xva, Yva = splittenths(X, Y, 8)
    #define parameter set for gradients
    p = params(model)
    #track the loss over epochs
    L = fill(NaN, epochs+1)
    L[1] = loss(Xva, Yva)
    for i ∈ 1:epochs
        Xtr, Ytr = shuffledata(Xtr, Ytr)
        #descend in minibatches
        for batch ∈ partition(1:size(Xtr,2), nbatch)
            grads = gradient(() -> loss(view(Xtr,:,batch), view(Ytr,:,batch)), p)
            update!(opt, p, grads)
        end
        L[i+1] = loss(Xva, Yva)
        println("$i) $(L[i+1])")
    end
    return L
end

## load and prepare CONUS dataframe

df = datadir("exp_pro", "conus_angles.csv") |> CSV.File |> DataFrame
#computes logslope, minorder, maxorder
derivedcols!(df)
#renames precip and temperature columns
renamePT!(df)

#columns of primary interest
predictors = [:angle, :logslope, :EVI, :maxorder, :minorder]
targets = [:P, :AI, :SSM, :SUSM]

#take a subset of the table
df = df[:,vcat(predictors, targets)]

## arrange the data for modeling

X = vcat(
    (df[!,predictors[1:3]] |> Matrix)' |> Matrix,
    onehotbatch(df.minorder, df.minorder |> unique |> sort),
    onehotbatch(df.maxorder, df.maxorder |> unique |> sort)
)
Y = (df[!,targets] |> Matrix)' |> Matrix
X = standardize(UnitRangeTransform, X, dims=2)
Y = standardize(UnitRangeTransform, Y, dims=2)
X, Y = shuffledata(X, Y)
NP = size(X,1)
NT = size(Y,1)
M = size(df,1)

#split into test/train
Xtr, Ytr, Xte, Yte = splittenths(X, Y, 8)

## first see how simple linear regression does

linmodel = Dense(NP, NT)
linloss(X, Y) = mse(linmodel(X), Y)

##

opt = Descent(0.1)
epochs = 6
nbatch = 5_000

path = train(Xtr, Ytr, linmodel, linloss, opt, epochs, nbatch)

println("linear test loss: $(linloss(Xte, Yte))")
plot(path)
xlabel("epoch")
ylabel("validation MSE loss");

## now see how a simple nonlinear network does

nnmodel = Chain(
    Dense(NP, 128, σ),
    Dense(128, 128, σ),
    Dense(128, NT)
)
nnloss(X, Y) = mse(nnmodel(X), Y)

##

opt = Descent(0.1)
epochs = 6
nbatch = 5_000

path = train(Xtr, Ytr, nnmodel, nnloss, opt, epochs, nbatch)

println("nn test loss: $(nnloss(Xte, Yte))")
plot(path)
xlabel("epoch")
ylabel("validation MSE loss");