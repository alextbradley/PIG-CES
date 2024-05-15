include("../shared.jl")


"""
output the truth in jld2 format using the data in the csv files

"""


#read in the data from the csv truth files
truth = CSV.read(joinpath(@__DIR__,"truth.csv"), DataFrame, header = false)
truth = Matrix(truth)';
noise = CSV.read(joinpath(@__DIR__,"noise.csv"), DataFrame, header = false)
noise = Matrix(noise)[1];
noise = map(Float64, noise);

Γ = noise * I #noisy observation of grounding line position
noise_dist = MvNormal(zeros(2), Γ) #two constraints
rng_seed = 12452
rng_model = Random.MersenneTwister(rng_seed)

y = truth .+ rand(rng_model, noise_dist)

@save "truth.jld2" y Γ rng_model