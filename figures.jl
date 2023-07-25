using Distributions, Random, Turing
using DataFrames, DataFramesMeta

prng = MersenneTwister(240819)
μ_true = 300; σ_true = 30; yield_true = Normal(μ_true, σ_true) |> x -> truncated(x, 0, Inf)

# Generate data
n_tests = 1

function get_test_data(n_tests::Int64, model = yield_true)
    prng = MersenneTwister(240819)
    return rand(prng, model, n_tests)
end

μ_pr = Normal(350, 50) |> x -> truncated(x, 0, Inf)
σ_pr = Exponential(10)

prng = MersenneTwister(240819)
prior_df = DataFrame(μ_pr = rand(prng, μ_pr, 10^3), σ = rand(prng, σ_pr, 10^3)) |>
    x -> @rtransform(x, :prior_samples = Normal(:μ_pr, :σ) |> 
                                            x -> truncated(x, 0, Inf) |> 
                                            x -> rand(prng, x, 1) |> 
                                            x -> x[1])

using Plots

histogram(prior_df.prior_samples, bins = 50, xlabel = "Yield strength (MPa)", 
          ylabel = "Frequency", legend = false, title = "Prior distribution")

# Define model
@model function yield_model(yield_strength_meas::Vector{Float64}, ϵ::Float64 = 1.0)
    
    # Priors
    σ ~ Exponential(10)
    μ ~ Normal(350, 50)
    
    # Gaussian model
    yield_strength ~ Normal(μ, σ)

    # Relating yield strength to imprecise test data
    n_samples = length(yield_strength)
    for n ∈ n_samples
        yield_strength_meas[n] ~ Normal(yield_strength, ϵ)
    end

end

# Set up MCMC sampling parameters

n_chains = 4; n_draws = 10^3; n_warmup = 10^3; thinning = 1; target_acceptance = 0.65
mcmc_sampler = Turing.NUTS(n_warmup, target_acceptance)

function posterior_samples(yield_data::Vector{Float64})
    return yield_model(yield_data) |>
        x -> sample(x, mcmc_sampler, MCMCThreads(),
                    n_draws, n_chains, progress = true, seed = 2408) |>
        x -> DataFrame(x) |>
        x -> @rtransform(x, :n_tests = length(yield_data)) |>
        x -> @select(x, :n_tests, :chain, :iteration, :μ, :σ)
end

yield_df = DataFrame(n_tests = Int64[], μ = Float64[], σ = Float64[])
for n in [1, 3, 5, 10, 50]
    yield_data = get_test_data(n)
    append!(yield_df, 
            posterior_samples(yield_data))
end

CSV.write("yield_df.csv", yield_df)

pwd()
