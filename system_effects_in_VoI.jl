#
# "System Effects in Identifying Risk-Optimal Data Requirements for Digital Twins of Structures"
#  Paper submitted to the journal of Reliability Engineering & System Safety
#  Domenic Di Francesco, PhD, CEng (MIMechE)
#  The Alan Turing Institute, University of Cambridge
#  
#  Accompanying Julia code to set-up and run calculations
#

######################################################
#
# Loading libraries
#
######################################################

# For describing probabilistic models
using Distributions, Turing, Random, LatinHypercubeSampling, Copulas
# For describing and solving decision problem
using JuMP, HiGHS, DecisionProgramming, LinearAlgebra
# For working with data
using CSV, DataFrames, DataFramesMeta


# For describing probabilistic models
using Distributions, Turing, Random, LatinHypercubeSampling, Copulas
# For describing and solving decision problem
using JuMP, HiGHS, DecisionProgramming, LinearAlgebra
# For working with data
using CSV, DataFrames, DataFramesMeta

######################################################
#
# Required functions
#
######################################################

function D_7608(S_MPa::Float64)
    return (Normal(log(10, 3.988 * 10^12) - 3*log(10, S_MPa), 0.2095))
end

function F2_7608(S_MPa::Float64)
    return (Normal(log(10, 1.231 * 10^12) - 3*log(10, S_MPa), 0.2279))
end

function draw_lhs(dist, n::Int; reprod::Int = 240819)
    Random.seed!(reprod)
    samples = randomLHC(n + 2, 1) |>
        x -> scaleLHC(x, [(0, 1)]) |>
        x -> quantile(dist, x)[:,1] |>
        x -> filter(!∈((-Inf, Inf)), x) |>
        x -> [x[i] for i in 1:length(x) if abs(x[i]) >= 10^-10]
    return samples
end

function get_Norm_param(log_μ::Float64, log_σ::Float64)
    mean = exp(log_μ + 1/2 * log_σ^2)
    sd = exp(log_μ + 1/2 * log_σ^2) * sqrt(exp(log_σ^2) - 1)
    return (μ = mean, σ = sd)
end

function get_LogNorm_params(μ::Float64, σ::Float64)
    log_sd = √(log(1 + σ^2 / μ^2))
    log_mean = log(μ) - 1/2 * log_sd^2
    return (log_μ = log_mean, log_σ = log_sd)
end

function cov_mat_2d(ρ::Float64, σ₁::Float64, σ₂::Float64)
    return([σ₁^2 ρ*σ₁*σ₂; ρ*σ₁*σ₂ σ₂^2])
end

n_baseline = (5 * 52 * 18 * 10) + (2 * 52 * 18 * 8)

function get_PoF(;n_years::Int = 1, n_cycles::Int = n_baseline,
                stress_samples::Vector{Float64} = prior_samples_df.stress,
                strength_samples::Vector{Float64} = prior_samples_df.σY,
                SCF_samples::Vector{Float64} = prior_samples_df.SCF,
                strengthen::Bool = false, βₛ::Float64 = 1.0, 
                replace::Bool = false, γᵣ::Float64 = 1.0, 
                reduce_op::Bool = false, βᵣ::Float64 = 1.0)
    
    if (replace == true)
        γᵣ = γᵣ * 0
    end

    SCF = SCF_samples .- (SCF_samples .- 1) * (1 - γᵣ)
    x_pred = stress_samples .* SCF; sort!(x_pred)

    if (strengthen == true)
        βₛ = βₛ * 3/4
    end

    if (reduce_op == true)
        βᵣ = βᵣ * 1/2
    end

    μ_pred = SN_post_df.C .- SN_post_df.m .* log.(10, x_pred * βₛ); y_pred = []
    for i in 1:n_samples
        append!(y_pred, 
                draw_lhs(Normal(μ_pred[i], SN_post_df.σ[i]), 100) |> x -> reduce(vcat, x))
    end

    PoF_SN = count(y_pred .<= log(10, n_years * n_cycles * βᵣ)) / length(y_pred)
    PoF_Ex = count(x_pred * βₛ .>= strength_samples) / length(x_pred)

    ΔPoF = PoF_SN + PoF_Ex - (PoF_SN * PoF_Ex)

    return(DataFrame(ΔPoF = minimum([1, (1 - (1 - ΔPoF)^n_years)]), βₛ = βₛ, βᵣ = βᵣ, γᵣ = γᵣ))
end

function gen_PoF_Dict(n; stress_samples::Vector{Float64}=prior_samples_df.stress, strength_samples::Vector{Float64}=prior_samples_df.σY, SCF_samples::Vector{Float64}=prior_samples_df.SCF)
    PoF_dict = Dict()
    args = Dict(:n_years => n, :stress_samples => stress_samples, 
                :strength_samples => strength_samples, :SCF_samples => SCF_samples)

    Threads.@threads for (k,v) in (("no_action", Dict()), ("strengthen", Dict(:strengthen => true)), 
                  ("replace", Dict(:replace => true)), ("reduce_operation", Dict(:reduce_op => true)), 
                  ("strengthen_replace", Dict(:strengthen => true, :replace => true)), 
                  ("strengthen_reduce", Dict(:strengthen => true, :reduce_op => true)), 
                  ("replace_reduce", Dict(:replace => true, :reduce_op => true)), 
                  ("strengthen_replace_reduce", Dict(:strengthen => true, :replace => true, :reduce_op => true)))
        v = merge(args, v)
        PoF_dict[k] = get_PoF(; v...).ΔPoF[1]
    end

    return PoF_dict
end

######################################################
#
# Bayesian models
#
######################################################

@model function SN_model(; N_meas::Int, log_S::Vector{Float64}, log_N::Vector{Float64})

    # Priors
    C ~ Normal(13, 3)
    m ~ truncated(Normal(3, 2), lower = 0)
    σ ~ Exponential(1)

    # Likelihood
    for n in 1:N_meas
        log_N[n] ~ Normal(C - m * log_S[n], σ)
    end

end

######################################################
#
# Setting up inputs
#
######################################################

sim_SN_data_df = DataFrame(); prng = MersenneTwister(240819)
for S in 80:5:120
    S = convert(Float64, S)
    append!(sim_SN_data_df, 
            DataFrame(S_MPa = S, log_N = D_7608(S) |> x -> rand(prng, x, 1)))
end

n_samples = 10_000; n_chains = 4; n_mcmc = Int(n_samples/n_chains)

SN_post_df = SN_model(N_meas = nrow(sim_SN_data_df), 
                      log_S = log.(10, sim_SN_data_df.S_MPa), 
                      log_N = sim_SN_data_df.log_N) |>
    x -> sample(MersenneTwister(240819), x, NUTS(), MCMCThreads(), n_mcmc, n_chains) |>
    x -> DataFrame(x) |>
    x -> @select(x, :iteration, :chain, :C, :m, :σ)

μ_l = Normal(50, 5); σ_l_params = get_LogNorm_params(10.0, 3.0)
μ_str = Normal(400, 20); σ_str_params = get_LogNorm_params(10.0, 3.0)
α_scf = Normal(2, 1/2) |> x -> truncated(x, lower = 0); θ_scf = Normal(1, 1/2) |> x -> truncated(x, lower = 0)

ρ_cop = 3/4

prior_df = DataFrame(
    μ_l = draw_lhs(μ_l, 100),
    σ_l = draw_lhs(LogNormal(σ_l_params.log_μ, σ_l_params.log_σ), 100),
    μ_str = draw_lhs(μ_str, 100),
    σ_str = draw_lhs(LogNormal(σ_str_params.log_μ, σ_str_params.log_σ), 100), 
    α_SCF = draw_lhs(α_scf, 100),
    θ_SCF = draw_lhs(θ_scf, 100)
) |>
    x -> @rtransform(x, :load = Normal(:μ_l, :σ_l)) |>
    x -> @rtransform(x, :load_samples = draw_lhs(:load, 100)) |>
    x -> @rtransform(x, :SCF_var = :α_SCF * :θ_SCF^2) |>
    x -> @rtransform(x, :G_cop = GaussianCopula(cov_mat_2d(ρ_cop, :σ_str, √:SCF_var))) |>
    x -> @rtransform(x, :SD = SklarDist(:G_cop, 
                                        (truncated(Normal(:μ_str, :σ_str), lower = 0), 
                                        censored(Gamma(:α_SCF, :θ_SCF), lower = 1)))) |>
    x -> @rtransform(x, :strength_samples = rand(MersenneTwister(240819), :SD, 100)[1, :]) |>
    x -> @rtransform(x, :SCF_samples = rand(MersenneTwister(240819), :SD, 100)[2,:])

prior_samples_df = DataFrame(
    stress = reduce(vcat, prior_df.load_samples),
    σY = reduce(vcat, prior_df.strength_samples),
    SCF = reduce(vcat, prior_df.SCF_samples)
)

n_years = 3; site_visit = 0.005; cost_strengthen = 0.015; cost_replace = 0.05; cost_reduce = 0.05

maint_opts = Dict("no_action" => 0, "strengthen" => cost_strengthen + site_visit, "replace" => cost_replace + site_visit, 
                  "reduce_operation" => cost_reduce, 
                  "strengthen_replace" => cost_strengthen + cost_replace + site_visit, 
                  "strengthen_reduce" => cost_strengthen + cost_reduce + site_visit, 
                  "replace_reduce" => cost_replace + cost_reduce + site_visit,
                  "strengthen_replace_reduce" => cost_strengthen + cost_replace + cost_reduce + site_visit)

maint_states = keys(maint_opts) |> x -> collect(x)
maint_values = [maint_opts[state] for state in maint_states]

β_states = ["Fail", "Survive"]; CoFs = [1, 0]

# strengthening reduces applied stress (e.g. by increasing area)
# replacing the component reduces the SCF due to misalignment to 1
# reducing operation reduces the number of stress cycles that the component sees

function solve_multi_year_id(;
    n_years::Int = 3,
    stress_samples::Vector{Float64} = prior_samples_df.stress,
    strength_samples::Vector{Float64} = prior_samples_df.σY,
    SCF_samples::Vector{Float64} = prior_samples_df.SCF,
    multi_thread = false)

    PoF_dict_y1 = gen_PoF_Dict(1,
                               stress_samples = stress_samples, 
                               strength_samples = strength_samples,
                               SCF_samples = SCF_samples)

    multi_year_decision = InfluenceDiagram()

    for i in 1:n_years
        add_node!(multi_year_decision, DecisionNode("maint$i", [], maint_states))
        add_node!(multi_year_decision, ValueNode("C_maint$i", ["maint$i"]))

        if i == 1
            add_node!(multi_year_decision, ChanceNode("β$i", ["maint$i"], β_states))
            add_node!(multi_year_decision, ValueNode("CoF$i", ["β$i"]))
        else
            deps = ["maint" * string(j) for j in 1:i]

            add_node!(multi_year_decision, ChanceNode("β$i", deps, β_states))
            add_node!(multi_year_decision, ValueNode("CoF$i", ["β$i"]))
        end
    end

    generate_arcs!(multi_year_decision)

    for n in 1:n_years

        Cₘ = UtilityMatrix(multi_year_decision, "C_maint$n")
        Cₘ["no_action"] = maint_opts["no_action"]
        Cₘ["strengthen"] = maint_opts["strengthen"]
        Cₘ["replace"] = maint_opts["replace"]
        Cₘ["reduce_operation"] = maint_opts["reduce_operation"]
        Cₘ["strengthen_replace"] = maint_opts["strengthen_replace"]
        Cₘ["strengthen_reduce"] = maint_opts["strengthen_reduce"]
        Cₘ["replace_reduce"] = maint_opts["replace_reduce"]
        Cₘ["strengthen_replace_reduce"] = maint_opts["strengthen_replace_reduce"]

        Cᵣ = UtilityMatrix(multi_year_decision, "CoF$n")
        Cᵣ["Fail"] = 1
        Cᵣ["Survive"] = 0

        add_utilities!(multi_year_decision, "C_maint$n", Cₘ)
        add_utilities!(multi_year_decision, "CoF$n", Cᵣ)

    end

    β₁ = ProbabilityMatrix(multi_year_decision, "β1")

    for i in maint_states
        β₁[i, :] = [PoF_dict_y1[i] (1 - PoF_dict_y1[i])]
    end

    β₂ = ProbabilityMatrix(multi_year_decision, "β2")

    for i in maint_states, j in maint_states
        β₂[i, j, :] = [1-(1-PoF_dict_y1[i])*(1-PoF_dict_y1[j]) (1-PoF_dict_y1[i])*(1-PoF_dict_y1[j])]
    end

    β₃ = ProbabilityMatrix(multi_year_decision, "β3")

    for i in maint_states, j in maint_states, k in maint_states
        β₃[i, j, k, :] = [1-(1-PoF_dict_y1[i])*(1-PoF_dict_y1[j])*(1-PoF_dict_y1[k]) 
                            (1-PoF_dict_y1[i])*(1-PoF_dict_y1[j])*(1-PoF_dict_y1[k])]
    end

    add_probabilities!(multi_year_decision, "β1", β₁)
    add_probabilities!(multi_year_decision, "β2", β₂)
    add_probabilities!(multi_year_decision, "β3", β₃)

    generate_diagram!(multi_year_decision)

    # Define and run solver
    if multi_thread == true
        multi_year_model = Model(optimizer_with_attributes(HiGHS.Optimizer, "threads" => Threads.nthreads()))
        set_optimizer_attribute(multi_year_model, "threads", Threads.nthreads())
    else
        multi_year_model = Model()
        set_optimizer(multi_year_model, HiGHS.Optimizer)
    end
    set_silent(multi_year_model)

    z = DecisionVariables(multi_year_model, multi_year_decision)
    EC = expected_value(multi_year_model, multi_year_decision, 
                        PathCompatibilityVariables(multi_year_model, multi_year_decision, z, probability_cut = false))

    @objective(multi_year_model, Min, EC)
    optimize!(multi_year_model)

    # Process results
    Z = DecisionStrategy(z)
    U_dist = UtilityDistribution(multi_year_decision, Z)

    # print_decision_strategy(multi_year_decision, Z,  StateProbabilities(multi_year_decision, Z))
    
    opt_df = DataFrame()

    for n in 1:n_years
        colname = "a_opt$n"
        opt_df[!, colname] = [maint_states[argmax(Z.Z_d[n])]]
    end

    opt_df.u_opt = [LinearAlgebra.dot(U_dist.p, U_dist.u)]

    return(opt_df)

end