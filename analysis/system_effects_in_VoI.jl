#=
    "System Effects in Identifying Risk-Optimal Data Requirements for Digital Twins of Structures"
    Paper submitted to the journal of Reliability Engineering & System Safety
    
    Accompanying Julia code to set-up and run calculations
   
    Domenic Di Francesco, PhD, CEng (MIMechE)
    The Alan Turing Institute, University of Cambridge
=#

#= ***************************************************************************
Loading libraries
*************************************************************************** =#

# For describing probabilistic models
using Distributions, Turing, Random, LatinHypercubeSampling, Copulas
# For describing and solving decision problem
using JuMP, HiGHS, Gurobi, DecisionProgramming, LinearAlgebra
# For working with data
using CSV, DataFrames, DataFramesMeta

#= ***************************************************************************
Required functions
*************************************************************************** =#

# Find the number of allowable cycles at a specified (constant-amplitude) stress range
# using the SN curves from BS7608

function D_7608(S_MPa::Float64)
    return (Normal(log(10, 3.988 * 10^12) - 3*log(10, S_MPa), 0.2095))
end

function F2_7608(S_MPa::Float64)
    return (Normal(log(10, 1.231 * 10^12) - 3*log(10, S_MPa), 0.2279))
end

# Draw latin hypercube samples from a distribution

function draw_lhs(dist, n::Int; reprod::Int = 240819)
    Random.seed!(reprod)
    samples = randomLHC(n + 2, 1) |>
        lhc -> scaleLHC(lhc, [(0, 1)]) |>
        lhc -> quantile(dist, lhc)[:,1] |>
        q -> filter(!∈((-Inf, Inf)), q) |>
        q -> [q[i] for i ∈ 1:length(q) if abs(q[i]) >= 10^-10]
    return samples
end

# Find the equivalent normal distribution parameters from lognormal distribution parameters

function get_Norm_param(log_μ::Float64, log_σ::Float64)
    mean = exp(log_μ + 1/2 * log_σ^2)
    sd = exp(log_μ + 1/2 * log_σ^2) * sqrt(exp(log_σ^2) - 1)
    return (μ = mean, σ = sd)
end

# Find the equivalent lognormal distribution parameters from normal distribution parameters

function get_LogNorm_params(μ::Float64, σ::Float64)
    log_sd = √(log(1 + σ^2 / μ^2))
    log_mean = log(μ) - 1/2 * log_sd^2
    return (log_μ = log_mean, log_σ = log_sd)
end

# Generate an n-dimensional covariance matrix

function Σ_mat(n::Int64, ρᵩ::Vector{Float64}, σ::Vector{Float64})
    @assert length(ρᵩ) == n*(n-1)/2; @assert length(σ) == n

    ρ_mat = Matrix{Float64}(I, n, n); Σ = zeros(n, n); i = 1
    
    for col ∈ 1:n, row ∈ 1:n
        if row > col
            ρ_mat[row, col] = ρᵩ[i]
            i += 1
        elseif col > row
            ρ_mat[row, col] = ρ_mat[col, row]
        end
        Σ[row, col] = row == col ? ρ_mat[row, col] * σ[row]^2 : 
                                   ρ_mat[row, col] * σ[row] * σ[col]
    end

    return Σ
end

# Generate an 2-dimensional covariance matrix

function cov_mat_2d(ρ::Float64, σ₁::Float64, σ₂::Float64)
    return([σ₁^2 ρ*σ₁*σ₂; ρ*σ₁*σ₂ σ₂^2])
end

#= 
Number of stress cycles per year for a rail bridge, which operates:
  - 10 trains per hour, 18 hours per day, on weekdays
  - 8 trains per hour, 18 hours per day, on weekends
=#

n_baseline = (5 * 52 * 18 * 10) + (2 * 52 * 18 * 8)

#=

Find the probability of failure (PoF) in n years, conditional on available interventions

Effect of each mitigation:
 - strengthening reduces applied stress (e.g. by increasing area)
 - replacing the component reduces the SCF due to misalignment to 1
 - reducing operation reduces the number of stress cycles that the component sees

=#

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
    for i ∈ 1:n_samples
        append!(y_pred, 
                draw_lhs(Normal(μ_pred[i], SN_post_df.σ[i]), 100) |> x -> reduce(vcat, x))
    end

    PoF_SN = count(y_pred .<= log(10, n_years * n_cycles * βᵣ)) / length(y_pred)
    PoF_Ex = count(x_pred * βₛ .>= strength_samples) / length(x_pred)

    ΔPoF = PoF_SN + PoF_Ex - (PoF_SN * PoF_Ex)

    return(DataFrame(ΔPoF = minimum([1, (1 - (1 - ΔPoF)^n_years)]), βₛ = βₛ, βᵣ = βᵣ, γᵣ = γᵣ))
end

# Find the probability of failure (PoF) in n years, conditional on available interventions

function gen_PoF_Dict(n; stress_samples::Vector{Float64} = prior_samples_df.stress, 
                      strength_samples::Vector{Float64} = prior_samples_df.σY, 
                      SCF_samples::Vector{Float64}=prior_samples_df.SCF)

    PoF_dict = Dict()
    args = Dict(:n_years => n, :stress_samples => stress_samples, 
                :strength_samples => strength_samples, :SCF_samples => SCF_samples)

    Threads.@threads for (k,v) ∈ (("no_action", Dict()), ("strengthen", Dict(:strengthen => true)), 
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

#= ***************************************************************************
Bayesian models
*************************************************************************** =#

@model function SN_model(; N_meas::Int, log_S::Vector{Float64}, log_N::Vector{Float64})

    # Define the prior distributions
    C ~ Normal(13, 3)
    m ~ truncated(Normal(3, 2), lower = 0)
    σ ~ Exponential(1)

    # Define the likelihood
    for n ∈ 1:N_meas
        log_N[n] ~ Normal(C - m * log_S[n], σ)
    end

end

@model function SHM(stress_meas::Float64, ϵ::Float64)
    # Define the prior distributions
    μ_stress ~ Normal(μ_l.μ , μ_l.σ)
    σ_stress ~ LogNormal(σ_l_params.log_μ, σ_l_params.log_σ)
    stress ~ Normal(μ_stress, σ_stress)

    # Define the likelihood
    stress_meas ~ Normal(stress, ϵ)
end

#= ***************************************************************************
Setting up inputs
*************************************************************************** =#

#=
Simulate test data from BS7608 SN curves, 
...accounting for epistemic uncertainty, using a Bayesian model
=#

sim_SN_data_df = DataFrame(); prng = MersenneTwister(240819)
for S ∈ 80:5:120
    S = convert(Float64, S)
    append!(sim_SN_data_df, 
            DataFrame(S_MPa = S, log_N = D_7608(S) |> x -> rand(prng, x, 1)))
end

n_samples = 10_000; n_chains = 4; n_mcmc = Int(n_samples/n_chains)

SN_post_df = SN_model(N_meas = nrow(sim_SN_data_df), 
                      log_S = log.(10, sim_SN_data_df.S_MPa), 
                      log_N = sim_SN_data_df.log_N) |>
    model -> sample(MersenneTwister(240819), model, NUTS(), MCMCThreads(), n_mcmc, n_chains) |>
    posterior -> DataFrame(posterior) |>
    post_df -> @select(post_df, :iteration, :chain, :C, :m, :σ)

# Specify prior models for load, SCF, and yield strength

μ_l = Normal(50, 5); σ_l_params = get_LogNorm_params(6.0, 3.0)
μ_str = Normal(400, 20); σ_str_params = get_LogNorm_params(10.0, 3.0)
α_scf = Normal(2, 1/2) |> x -> truncated(x, lower = 0); γ_scf = Normal(1/2, 1/2) |> x -> truncated(x, lower = 0)

ρ_cop = 2/3; n_lhc = 100

prior_df = DataFrame(μ_l = draw_lhs(μ_l, n_lhc),
                     σ_l = draw_lhs(LogNormal(σ_l_params.log_μ, σ_l_params.log_σ), n_lhc),
                     μ_str = draw_lhs(μ_str, n_lhc),
                     σ_str = draw_lhs(LogNormal(σ_str_params.log_μ, σ_str_params.log_σ), n_lhc),
                     α_SCF = draw_lhs(α_scf, n_lhc),
                     γ_SCF = draw_lhs(γ_scf, n_lhc)) |>
    df -> @rtransform(df, :load = Normal(:μ_l, :σ_l)) |>
    df -> @rtransform(df, :load_samples = draw_lhs(:load, n_lhc)) |>
    df -> @rtransform(df, :SCF_var = :α_SCF * :γ_SCF^2) |>
    df -> @rtransform(df, :G_cop = GaussianCopula(cov_mat_2d(ρ_cop, :σ_str, √:SCF_var))) |>
    df -> @rtransform(df, :SD = SklarDist(:G_cop, 
                                         (truncated(Normal(:μ_str, :σ_str), lower = 0), 
                                         censored(Gamma(:α_SCF, :γ_SCF), lower = 1)))) |>
    df -> @rtransform(df, :strength_samples = rand(MersenneTwister(240819), :SD, n_lhc)[1, :]) |>
    df -> @rtransform(df, :SCF_samples = rand(MersenneTwister(240819), :SD, n_lhc)[2,:])

prior_samples_df = prior_df |>
    df -> DataFrame(stress = reduce(vcat, df.load_samples),
                    σY = reduce(vcat, df.strength_samples),
                    SCF = reduce(vcat, df.SCF_samples))


# Define the costs of interventions for use in the decision problem

n_years = 3; site_visit = 0.01; cost_strengthen = 0.025; cost_replace = 0.075; cost_reduce = 0.05

maint_opts = Dict("no_action" => 0, "strengthen" => cost_strengthen + site_visit, "replace" => cost_replace + site_visit, 
                  "reduce_operation" => cost_reduce, 
                  "strengthen_replace" => cost_strengthen + cost_replace + site_visit, 
                  "strengthen_reduce" => cost_strengthen + cost_reduce + site_visit, 
                  "replace_reduce" => cost_replace + cost_reduce + site_visit,
                  "strengthen_replace_reduce" => cost_strengthen + cost_replace + cost_reduce + site_visit)

maint_states = keys(maint_opts) |> x -> collect(x)
maint_values = [maint_opts[state] for state ∈ maint_states]

β_states = ["Fail", "Survive"]

CoFs = [1, 0]

#=
Create a function to solve the decision problem as an influence diagram
returning a dataframe of the expected optimal utility and associated decision
=#

function solve_id(;
    stress_samples::Vector{Float64} = prior_samples_df.stress,
    strength_samples::Vector{Float64} = prior_samples_df.σY,
    SCF_samples::Vector{Float64} = prior_samples_df.SCF)

    # Find PoFs from prior models and measurement data
    PoF_dict_y1 = gen_PoF_Dict(1,
                               stress_samples = stress_samples, 
                               strength_samples = strength_samples,
                               SCF_samples = SCF_samples)

    # Initialise influence diagram and add nodes
    SIM = InfluenceDiagram()

    add_node!(SIM, DecisionNode("maint", [], maint_states))
    add_node!(SIM, ChanceNode("β", ["maint"], β_states))
    add_node!(SIM, ValueNode("CoF", ["β"]))
    add_node!(SIM, ValueNode("C_maint", ["maint"]))

    generate_arcs!(SIM)

    # Populate structural reliability node(s) with PoFs
    β = ProbabilityMatrix(SIM, "β")
    
    for i ∈ maint_states
        β[i, :] = [PoF_dict_y1[i] (1 - PoF_dict_y1[i])]
    end

    # Populate maintenance cost node(s)
    Cₘ = UtilityMatrix(SIM, "C_maint")
    for i ∈ maint_states
        Cₘ[i] = maint_opts[i]
    end

    # Populate failre cost node(s)
    Cᵣ = UtilityMatrix(SIM, "CoF")
    Cᵣ["Fail"] = 1
    Cᵣ["Survive"] = 0

    add_probabilities!(SIM, "β", β)
    add_utilities!(SIM, "C_maint", Cₘ)
    add_utilities!(SIM, "CoF", Cᵣ)

    generate_diagram!(SIM)

    # Define JuMP model with solver using all available threads
    SIM_model = Model()
    set_optimizer(SIM_model, Gurobi.Optimizer) # or HiGHS.Optimizer if no Gurobi license is available
    set_optimizer_attribute(SIM_model, "threads", Threads.nthreads())

    # Define decision variables and expected utility for optimisation
    z = DecisionVariables(SIM_model, SIM)
    EC = expected_value(SIM_model, SIM, 
                        PathCompatibilityVariables(SIM_model, SIM, z))

    @objective(SIM_model, Min, EC)
    optimize!(SIM_model)

    # Extract a* and u* from the solution
    Z = DecisionStrategy(z)
    U_dist = UtilityDistribution(SIM, Z)
    
    # Return results as a dataframe
    opt_df = DataFrame(a_opt = maint_states[argmax(Z.Z_d[1])],
                       u_opt = LinearAlgebra.dot(U_dist.p, U_dist.u))

    return(opt_df)

end

#=
Create a function to solve the sequential decision problem as an influence diagram
returning a dataframe of the expected optimal utility and associated decision
=#

function solve_multi_year_id(;
    n_years::Int = 3,
    stress_samples::Vector{Float64} = prior_samples_df.stress,
    strength_samples::Vector{Float64} = prior_samples_df.σY,
    SCF_samples::Vector{Float64} = prior_samples_df.SCF)

    # Find PoFs from prior models and measurement data
    PoF_dict_y1 = gen_PoF_Dict(1,
                               stress_samples = stress_samples, 
                               strength_samples = strength_samples,
                               SCF_samples = SCF_samples)

    # Initialise influence diagram and add nodes
    multi_year_decision = InfluenceDiagram()

    for i ∈ 1:n_years
        add_node!(multi_year_decision, DecisionNode("maint$i", [], maint_states))
        add_node!(multi_year_decision, ValueNode("C_maint$i", ["maint$i"]))

        if i == 1
            add_node!(multi_year_decision, ChanceNode("β$i", ["maint$i"], β_states))
            add_node!(multi_year_decision, ValueNode("CoF$i", ["β$i"]))
        else
            deps = ["maint" * string(j) for j ∈ 1:i]

            add_node!(multi_year_decision, ChanceNode("β$i", deps, β_states))
            add_node!(multi_year_decision, ValueNode("CoF$i", ["β$i"]))
        end
    end

    generate_arcs!(multi_year_decision)

    # Populate maintenance and failure cost node(s)
    for n ∈ 1:n_years

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

    # Populate structural reliability node with PoFs
    β₁ = ProbabilityMatrix(multi_year_decision, "β1")
    for i ∈ maint_states
        β₁[i, :] = [1-(1-PoF_dict_y1[i])
                      (1 - PoF_dict_y1[i])]
    end

    β₂ = ProbabilityMatrix(multi_year_decision, "β2")
    for i ∈ maint_states, j ∈ maint_states
        β₂[i, j, :] = [1-(1-PoF_dict_y1[i])*(1-PoF_dict_y1[j])
                         (1-PoF_dict_y1[i])*(1-PoF_dict_y1[j])]
    end

    β₃ = ProbabilityMatrix(multi_year_decision, "β3")
    for i ∈ maint_states, j ∈ maint_states, k ∈ maint_states
        β₃[i, j, k, :] = [1-(1-PoF_dict_y1[i])*(1-PoF_dict_y1[j])*(1-PoF_dict_y1[k]) 
                            (1-PoF_dict_y1[i])*(1-PoF_dict_y1[j])*(1-PoF_dict_y1[k])]
    end

    add_probabilities!(multi_year_decision, "β1", β₁)
    add_probabilities!(multi_year_decision, "β2", β₂)
    add_probabilities!(multi_year_decision, "β3", β₃)

    generate_diagram!(multi_year_decision)

    # Define JuMP model with solver using all available threads
    multi_year_model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(multi_year_model, "threads", Threads.nthreads())

    # Define decision variables and expected utility for optimisation
    z = DecisionVariables(multi_year_model, multi_year_decision)
    EC = expected_value(multi_year_model, multi_year_decision, 
                        PathCompatibilityVariables(multi_year_model, multi_year_decision, z))

    @objective(multi_year_model, Min, EC)
    optimize!(multi_year_model)

    # Extract a* and u* from the solution
    Z = DecisionStrategy(z)
    U_dist = UtilityDistribution(multi_year_decision, Z)

    # Return results as a dataframe
    opt_df = DataFrame()
    for n ∈ 1:n_years
        colname = "a_opt$n"
        opt_df[!, colname] = [maint_states[argmax(Z.Z_d[n])]]
    end

    opt_df.u_opt = [LinearAlgebra.dot(U_dist.p, U_dist.u)]

    return(opt_df)

end

#=

Running VoPI analysis

=#

function VoPI(;data::String, n_step::Int64 = 10)
    data_types = ["none", "inspection", "testing", "SHM", "inspection + testing", 
                  "inspection + SHM", "testing + SHM", "inspection + testing + SHM"]
    @assert data ∈ data_types

    iter_df = Vector{DataFrame}(undef, nrow(prior_samples_df) ÷ n_step)

    if data == "none"
        VoPI_df = solve_id() |>
            df -> @rtransform(df, :stress_meas = missing) |>
            df -> @rtransform(df, :strength_meas = missing) |>
            df -> @rtransform(df, :SCF_meas = missing)
    
    else
        progress = 0
        for i in 1:size(iter_df)
            inspection_data = occursin("inspection", data) ? [prior_samples_df.SCF[i]] : prior_samples_df.SCF
            test_data = occursin("testing", data) ? [prior_samples_df.σY[i]] : prior_samples_df.σY
            SHM_data = occursin("SHM", data) ? [prior_samples_df.stress[i]] : prior_samples_df.stress

            iter_df[i] = solve_id(SCF_samples = inspection_data, stress_samples = SHM_data, strength_samples = test_data) |>
                df -> @rtransform(df, :stress_meas = length(SHM_data) == 1 ? prior_samples_df.stress[i] : missing) |>
                df -> @rtransform(df, :strength_meas = length(test_data) == 1 ? prior_samples_df.σY[i] : missing) |>
                df -> @rtransform(df, :SCF_meas = length(inspection_data) == 1 ? prior_samples_df.SCF[i] : missing)
            progress += n_step
            println("\n\n", progress / nrow(prior_samples_df) * 100, "% complete\n\n")
        end

        VoPI_df = DataFrame()
        for df in iter_df
            append!(VoPI_df, df)
        end
    end

   return VoPI_df 
end

# prior decision analysis
prior_decision = VoPI(data = "none")

#=
Running sensitivity analysis, demonstrating effect of measurement uncertainty on VoI
=# 

n_samples = 10_000; n_chains = 4; n_mcmc = Int(n_samples/n_chains); VoSHM_df = DataFrame()

for meas_error ∈ [1.0, 5.0, 10.0, 20.0]
    progress = 1
    for i ∈ 1:length(VoSHM_iter_df)
        σ_meas = prior_samples_df.stress[i]
        SHM(σ_meas, meas_error) |>
            # Sample from joint posterior π(θ|z) of load model conditional on sensor data
            model -> sample(MersenneTwister(240819), model, NUTS(), 
                            MCMCThreads(), n_mcmc, n_chains) |>
            posterior -> DataFrame(posterior) |>
            #= draw samples from posterior predicitve distribution of stress
            π(zₚ|z) = ∫ π(zₚ|θ) × π(θ|z) dθ =#
            posterior_df -> reduce(vcat, posterior_df.stress) |>
            #= Propagate through influence diagram and find Exp. optimal utility (and intervention)
            a* = argmax_{a ∈ A} E[u(a, π(θ∣z))] =#
            σ -> solve_id(stress_samples = σ) |>
            df -> @rtransform(df, :meas = σ_meas, error = meas_error) |>
            df -> append!(VoSHM_df, df)
        progress += 1
        println("\n\n Running VoSHM with error: $meas_error \n", 
                100 * progress / length(VoSHM_iter_df) |> x -> round(x, digits = 2), "% complete")
    end
end