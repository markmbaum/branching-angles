using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using CSV
using DataFrames
using Flux
using Flux: mse, params, update!
using Base.Iterators: partition
using PyPlot
using Random: shuffle
using StatsBase
using Statistics

pygui(true)

## load and prepare CONUS dataframe

df = datadir("exp_pro", "conus_angles.csv") |> CSV.File |> DataFrame
#computes logslope, minorder, maxorder
derivedcols!(df)
#renames precip and temperature columns
renamePT!(df)
#columns of primary interest
cols = [
    :P,
    :T,
    :EVI,
    :SSM,
    :SUSM,
    :AI,
    :logslope,
    :maxorder,
    :minorder
]
#take a subset of the most interestig columns
df = df[:,vcat(:angle, cols)]
#standardize all the numeric quantities
transform!(df, cols.=>zscore.=>cols, :angle=>zscore=>:angle)

##

function shuffledata(X, y)
    @assert size(X,2) == size(y,2)
    idx = shuffle(1:size(X,2))
    X[:,idx], y[:,idx]
end

## arrange the data for modeling

N = length(cols)
M = size(df,1)
X = (df[!,cols] |> Matrix)' |> Matrix
y = (df.angle |> Vector)' |> Matrix
X, y = shuffledata(X, y)

## a data splitting function

function splittenths(X, y, tenth)
    M = size(X,2)
    idx = tenth*(M ÷ 10)
    X[:,1:idx], y[:,1:idx], X[:,idx+1:end], y[:,idx+1:end]
end

## a model training function

function descend!(X, y, model, loss, opt, epochs, nbatch)
    #data shape
    N, M = size(X)
    @assert size(y,2) == M
    #split for train/validation
    Xtr, ytr, Xva, yva = splittenths(X, y, 9)
    #define parameter set for gradients
    p = params(model)
    #track the loss over epochs
    L = fill(NaN, epochs+1)
    L[1] = loss(Xva, yva)
    for i ∈ 1:epochs
        Xtr, ytr = shuffledata(Xtr, ytr)
        #descend in minibatches
        for batch ∈ partition(1:size(Xtr,2), nbatch)
            grads = gradient(() -> loss(view(Xtr,:,batch), view(ytr,:,batch)), p)
            update!(opt, p, grads)
        end
        L[i+1] = loss(Xva, yva)
    end
    return L
end

## first see how simple linear regression does

model = Dense(N, 1)
loss(x, y) = mse(model(x), y)
opt = Descent(5e-2)
epochs = 10
nbatch = 10_000 #lots of batches go faster
Xtr, ytr, Xte, yte = splittenths(X, y, 8)

path = descend!(Xtr, ytr, model, loss, opt, 5, nbatch)

println("test loss: $(loss(Xte, yte))")
plot(path)
xlabel("epoch")
ylabel("validation MSE loss")

## now see how a simple nonlinear network does

model = Chain(
    Dense(N, 2N, σ),
    Dense(2N, 2N, σ),
    Dense(2N, 1)
)
loss(x, y) = mse(model(x), y)
opt = Descent(5e-2)
epochs = 10
nbatch = 10_000 #lots of batches go faster
Xtr, ytr, Xte, yte = splittenths(X, y, 8)

path = descend!(Xtr, ytr, model, loss, opt, 5, nbatch)

println("test loss: $(loss(Xte, yte))")
plot(path)
xlabel("epoch")
ylabel("validation MSE loss")