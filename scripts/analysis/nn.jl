using DrWatson
@quickactivate "Branching Angles"
push!(LOAD_PATH, srcdir())
using BranchingAngles
using Arrow
using DataFrames
using Flux
using Random: shuffle
using ProgressMeter: @showprogress
using PyPlot
using StatsBase: mean, zscore

pygui(true)

##

#function to split a matrix into training, validation, and test sets *randomly*
function split(X, Y, validation::Int=1, test::Int=1)
    n = size(X,2)
    @assert size(Y,2) == n
    idx = shuffle(1:n)
    t = n ÷ 10
    tr = t*(10 - validation - test)
    va = t*(10 - validation)
    return(
        X[:,idx[1:tr]], #training data
        Y[:,idx[1:tr]], #training labels
        X[:,idx[tr+1:va]], #validation data
        Y[:,idx[tr+1:va]], #validation labels
        X[:,idx[va+1:end]], #test data
        Y[:,idx[va+1:end]], #test labels 
    )
end

function train(model, loss, Xtr, Ytr, Xva, Yva, Xte, Yte; nepoch=8, nbatch=32, verbose=false)
    #define training strategy
    dl = Flux.DataLoader((data=Xtr, labels=Ytr), batchsize=nbatch)
    p = Flux.params(model)
    opt = ADAM()
    #train the model
    history = (tr_loss=zeros(nepoch+1), va_loss=zeros(nepoch+1))
    history[:tr_loss][1] = loss(Xtr, Ytr)
    history[:va_loss][1] = loss(Xva, Yva)
    @showprogress for i ∈ 1:nepoch
        if verbose
            println("epoch $i")
            println("  train loss: $(history[:tr_loss][i])")
            println("  valid loss: $(history[:va_loss][i])")
        end
        Flux.train!(loss, p, dl, opt)
        history[:tr_loss][i+1] = loss(Xtr, Ytr)
        history[:va_loss][i+1] = loss(Xva, Yva)
    end
    te_loss = loss(Xte, Yte)
    verbose && println("test loss: $te_loss")
    return(history, te_loss)
end

## load and prepare CONUS data

df = datadir("exp_pro", "conus_angles.feather") |> Arrow.Table |> DataFrame
#computes logslope, minorder, maxorder
derivedcols!(df)
#renames precip and temperature columns
renamePT!(df)

## bin by watershed and take the means

inputs = [:angle, :logslope, :minorder, :maxorder]
labels = [:P, :T, :AI, :EVI, :SSM]

ge = combine(
    groupby(
        df,
        :huc8
    ),
    inputs .=> mean .=> inputs,
    labels .=> mean .=> labels,
    nrow => :n
)
filter!(r -> r.n > 10, ge)
transform!(
    ge,
    inputs .=> zscore .=> inputs,
    labels .=> zscore .=> labels
)

##

arrdata(df) = df |> Matrix |> transpose |> Matrix .|> Float32

X = ge[!,inputs] |> arrdata
Y = ge[!,labels] |> arrdata

Xtr, Ytr, Xva, Yva, Xte, Yte = split(X, Y)

##

np = 128
model = Chain(
    Dense(size(X,1) => np, relu),
    Dense(np => np, relu),
    Dense(np => np, relu),
    Dense(np => size(Y,1))
)

norm(x) = sum(abs2, x)
penalty() = sum(norm, Flux.params(model))
loss(X, Y) = Flux.mae(model(X), Y) + 0.01*penalty()
loss(Z) = loss(Z...)

##

history, te_loss = train(model, loss, Xtr, Ytr, Xva, Yva, Xte, Yte, nepoch=48)