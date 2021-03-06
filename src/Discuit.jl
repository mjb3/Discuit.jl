"""

Discuit is a package for:

- Customisable finite adaptive MCMC algorithm for fast parameter inference.
- Pooley model based proposal (MBP) method for improved mixing.
- Simulation via the Gillespie direct method.
- Automated tools for convergence diagnosis and analysis.
- Developed for Julia `1.0`.
- Author: Martin Burke (martin.burke@bioss.ac.uk)
- Date: 2018-08-22

"""
module Discuit

## depends CHANGE TO IMPORT
import Distributions
import Random
import CSV
import DataFrames

## exports
# public structs
export DiscuitModel, Trajectory, Observations, SimResults, MCMCResults, MultiMCMCResults, AutocorrelationResults
# core functionality
export set_random_seed, gillespie_sim, run_single_chain_analysis, run_multi_chain_analysis, compute_autocorrelation
# model helpers
export generate_model, generate_gaussian_obs_model, generate_generic_obs_function, generate_weak_prior
# utilities (e.g. saving results to file0
export print_trajectory, print_observations, print_results, print_autocorrelation, get_observations, tabulate_results
# visualisation
export plot_trajectory, plot_parameter_trace, plot_parameter_marginal, plot_parameter_heatmap, plot_geweke_series, plot_autocorrelation
# custom functionality (in development)
export MarkovState, ParameterProposal, PrivateDiscuitModel, generate_custom_x0, run_custom_single_chain_analysis, run_custom_multi_chain_analysis
# for unit testing:
# export PrivateDiscuitModel, get_private_model, model_based_proposal, standard_proposal, gillespie_sim_x0, run_geweke_test, compute_full_log_like # Event


include("./discuit_structs.jl")
include("./discuit_models.jl")
include("./discuit_visualisation_uc.jl")

## constants
# max trajectory constant (make this optional?)
const MAX_TRAJ = 196000
const DF_PROP_LIKE = 1.0
const NULL_LOG_LIKE = -Inf
const AC_LAG_INT = 10       # number of autocorrelation lag intervals
const C_DEBUG = true

## model generating functions
# TO BE ADDED *******

## FOR DEBUG
"""
    set_random_seed(seed)

# Examples

    set_random_seed(1234)

Does what it says on the tin but only if you give it an integer.
"""
function set_random_seed(seed::Int64)
    Random.seed!(seed);
end

## gillespie sim
# choose event type
function choose_event(cum_rates::Array{Float64,1})
    etc = rand() * cum_rates[end]
    for i in 1:(length(cum_rates) - 1)
        cum_rates[i] > etc && return i
    end
    return length(cum_rates)
end
# iterate particle
function iterate_sim!(model::PrivateDiscuitModel, trajectory::Trajectory, population::Array{Int64}, parameters::Array{Float64,1}, time::Float64, tmax::Float64)
    # declare array for use by loop
    cum_rates = Array{Float64, 1}(undef, size(model.m_transition, 1))
    while true
        model.rate_function(cum_rates, parameters, population)
        cumsum!(cum_rates, cum_rates)
        # get_cum_rates!(cum_rates, event_types, parameters, population)
        # 0 rate test
        cum_rates[end] == 0.0 && break
        time -= log(rand()) / cum_rates[end]
        # break if max time exceeded
        time > tmax && break
        # else choose event type (init as final event)
        et = choose_event(cum_rates)
        # add event to trajectory
        push!(trajectory.time, time)
        push!(trajectory.event_type, et)
        # update population
        population .+= model.m_transition[et,:] # event_types[et].v_transition
        # println(" t: ", time, ". R: ", cum_rates[end])
    end
end

## run sim and return trajectory
# - ADD option for >1 sim?
"""
    gillespie_sim(model, parameters, tmax = 100.0, num_obs = 5)

**Parameters**
- `model`       -- `DiscuitModel` (see [Discuit.jl models]@ref).
- `parameters`  -- model parameters.
- `tmax`        -- maximum time.
- `num_obs`     -- number of observations to draw,

Run a Doob-Gillespie simulation on `model`. Returns a SimResults containing the trajectory and observations data.
"""
function gillespie_sim(model::DiscuitModel, parameters::Array{Float64,1}, tmax::Float64 = 100.0, num_obs::Int64 = 5)
    # initialise some things
    p_model = get_private_model(model, Observations([0.0], [0 0]) ) # NOTE TO JAMIE: can this be made nullable? or is it better to use delegates?
    obs_times = collect(tmax / num_obs : tmax / num_obs : tmax)
    obs_vals = Array{Int64, ndims(p_model.initial_condition) + 1}(undef, length(obs_times), length(p_model.initial_condition))
    # time = 0.0
    population = copy(p_model.initial_condition)
    trajectory = Trajectory(Float64[], Int64[])
    # run
    println("running DGA simulation...")
    t_prev = model.t0_index == 0 ? 0.0 : parameters[model.t0_index]
    for i in eachindex(obs_times)
        iterate_sim!(p_model, trajectory, population, parameters, t_prev, obs_times[i])
        obs_vals[i,:] .= p_model.observation_function(population)
        t_prev = obs_times[i]
    end
    println(" finished (", length(trajectory.time), " events).")
    # print trajectory
    population .= p_model.initial_condition
    pop_long = zeros(Int64, length(trajectory.time), length(population))
    for i in eachindex(trajectory.time)
        population .+= p_model.m_transition[trajectory.event_type[i], :]
        pop_long[i,:] .= population
    end
    # return trajectory
    return SimResults(model.model_name, trajectory, pop_long, Observations(obs_times, obs_vals))
end
# sim to initialise Markov chain
function gillespie_sim_x0(model::PrivateDiscuitModel, parameters::Array{Float64,1}, full_like::Bool)
    # generate = true
    while true
        # initialise some things
        population = copy(model.initial_condition)
        trajectory = Trajectory(Float64[], Int64[])
        output = 0.0    # obs model
        # run
        t_prev = model.t0_index == 0 ? 0.0 : parameters[model.t0_index]
        # println(" t0: ", t_prev)
        for i in eachindex(model.obs_data.time)
            # println(" period: ", i, ". time: ", t_prev)
            iterate_sim!(model, trajectory, population, parameters, t_prev, model.obs_data.time[i])
            # evaluate log likelihood
            output += model.observation_model(model.obs_data.val[i,:], population)
            output == NULL_LOG_LIKE && break
            t_prev = model.obs_data.time[i]
        end
        ## REPLACE WITH CONST OR PARAMETER? *****
        if output != NULL_LOG_LIKE # && length(trajectory.time) > e.g. 5?
            # return MarkovState(ParameterProposal(parameters, model.prior_density(parameters)), trajectory, full_like ? compute_full_log_like(model, parameters, trajectory) : output, DF_PROP_LIKE, 0)
            return MarkovState(ParameterProposal(parameters), trajectory, full_like ? compute_full_log_like(model, parameters, trajectory) : output, DF_PROP_LIKE, 0)
        end
    end
end

## generate MarkovState based on custom trajectory (helper)
"""
    generate_custom_x0(model, obs_data, parameters, event_times, event_types)

**Parameters**
- `model`               -- `DiscuitModel` (see [Discuit.jl models]@ref).
- `obs_data`            -- `Observations` data.
- `parameters`          -- `Array` of initial model parameters.
- `event_times`         -- `Array` of floats representing the event times.
- `event_types`         -- `Array` of integer event types.

Generate an initial `MarkovState` for use in a custom MCMC algorithm.
"""
function generate_custom_x0(model::DiscuitModel, obs_data::Observations, parameters::Array{Float64, 1}, event_times::Array{Float64, 1}, event_types::Array{Int, 1})
    prop = ParameterProposal(parameters)
    trajectory = Trajectory(event_times, event_types)
    ## return as Proposal result
    return MarkovState(prop, trajectory, compute_full_log_like(get_private_model(model, obs_data), parameters, trajectory), DF_PROP_LIKE, 0)
end

## mv parameter proposal
function get_mv_param(model::PrivateDiscuitModel, g::Distributions.MvNormal, sclr::Float64, theta_i::Array{Float64, 1})
    output = rand(g)
    output .*= sclr
    output .+= theta_i
    return ParameterProposal(output)
end

## standard trajectory proposal
# likelihood function
function compute_full_log_like(model::PrivateDiscuitModel, parameters::Array{Float64,1}, trajectory::Trajectory)
    population = copy(model.initial_condition)
    # log_like = 0.0
    ll_traj = 0.0
    ll_obs = 0.0
    t = 0.0
    if model.t0_index > 0
        trajectory.time[1] < parameters[model.t0_index] && (return -Inf)
        t = parameters[model.t0_index]
    end
    evt_i = 1
    # workspace
    lambda = Array{Float64, 1}(undef, size(model.m_transition, 1))
    # for each observation period:
    for obs_i in eachindex(model.obs_data.time)
        while evt_i <= length(trajectory.time)
            trajectory.time[evt_i] > model.obs_data.time[obs_i] && break
            model.rate_function(lambda, parameters, population)
            ## deal with event
            # event log likelihood
            ll_traj += log(lambda[trajectory.event_type[evt_i]]) - (sum(lambda) * (trajectory.time[evt_i] - t))
            # evt_i < 4 && println(" i: ", evt_i, ". ll = ln ", lambda[trajectory.event_type[evt_i]], " - ", sum(lambda), "*", trajectory.time[evt_i] - t ,"=", ll_traj)
            # update population and handle -ve (NB. template? check not req'd for MBP)
            population .+= model.m_transition[trajectory.event_type[evt_i],:]
            any(x->x<0, population) && return -Inf
            t = trajectory.time[evt_i]
            evt_i += 1
        end
        ## handle observation
        model.rate_function(lambda, parameters, population)
        # log likelihood
        ll_traj -= sum(lambda) * (model.obs_data.time[obs_i] - t)
        # obs model
        ll_obs += model.observation_model(model.obs_data.val[obs_i,:], population)
        ll_traj == NULL_LOG_LIKE && return log(0.0)
        t = model.obs_data.time[obs_i]
    end
    # model.fc_index > 0 && (population[model.fc_index] == 0 || (log_like = NULL_LOG_LIKE))
    # println(" ll: ", ll_traj, " ::: ", ll_obs)
    return ll_traj + ll_obs
end

## insert event at appropriate index (std)
function add_event!(xf_trajectory::Trajectory, evt_tp::Int64, evt_tm::Float64)
    if (length(xf_trajectory.time) == 0 || evt_tm > xf_trajectory.time[end])
        push!(xf_trajectory.time, evt_tm)
        push!(xf_trajectory.event_type, evt_tp)
    else
        for i in eachindex(xf_trajectory.time)
            if xf_trajectory.time[i] > evt_tm
                insert!(xf_trajectory.time, i, evt_tm)
                insert!(xf_trajectory.event_type, i, evt_tp)
                break
            end
        end
    end
end

## event type count
function get_event_type_count(trajectory::Trajectory, et::Int64)
    output = 0
    for i in eachindex(trajectory.event_type)
        trajectory.event_type[i] == et && (output += 1)
    end
    return output
end

## standard proposal function (VALIDATED)
function standard_proposal_old(model::PrivateDiscuitModel, xi::MarkovState, xf_parameters::ParameterProposal)
    ## choose proposal type
    prop_type = rand(1:3)
    # trajectory proposal
    # - NEED TO MAKE THIS MORE EFFICIENT ****
    xf_trajectory = deepcopy(xi.trajectory)
    t0 = (model.t0_index == 0) ? 0.0 : xf_parameters.value[model.t0_index]
    if prop_type == 3
        ## move
        length(xi.trajectory.time) == 0 && (return MarkovState(xf_parameters, xi.trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, prop_type))
        # - IS THERE A MORE EFFICIENT WAY TO DO THIS? I.E. ROTATE using circshift or something?
        # choose event and define new one
        evt_i = rand(1:length(xi.trajectory.time))
        # evt_tm = (rand() * (model.obs_data.time[end] - t0)) + t0 #, xi.trajectory.event_type[evt_i])
        evt_tp = xi.trajectory.event_type[evt_i]
        # remove old one
        splice!(xf_trajectory.time, evt_i)
        splice!(xf_trajectory.event_type, evt_i)
        # add new one at random time
        add_event!(xf_trajectory, evt_tp, (rand() * (model.obs_data.time[end] - t0)) + t0)
        # compute ln g(x)
        prop_lk = 1.0
    else
        ## choose type and count
        tp = rand(1:size(model.m_transition, 1))
        ec = get_event_type_count(xi.trajectory, tp)
        ## insert / delete
        if prop_type == 1
            ## insert at randomly chosen time
            add_event!(xf_trajectory, tp, (rand() * (model.obs_data.time[end] - t0)) + t0)
            ## compute ln g(x)
            prop_lk = (model.obs_data.time[end] - t0) / (ec + 1)
        else
            ## delete
            # println(" deleting... tp:", tp, " - ec: ", ec)
            ec == 0 && (return MarkovState(xi.parameters, xf_trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, prop_type))
            # choose event index (CAN MAKE THIS MORE EFFICIENT ***)
            evt_t_i = rand(1:ec)
            evt_t_c = 0
            for i in eachindex(xf_trajectory.event_type)
                xf_trajectory.event_type[i] == tp && (evt_t_c += 1)
                if evt_t_i == evt_t_c
                    # remove
                    splice!(xf_trajectory.time, i)
                    splice!(xf_trajectory.event_type, i)
                    break
                end
            end
            # remove
            # splice!(xf_trajectory.time, evt_i)
            # splice!(xf_trajectory.event_type, evt_i)
            # compute ln g(x)
            prop_lk = ec / (model.obs_data.time[end] - t0)
        end # end of insert/delete
    end
    ## evaluate full likelihood for trajectory proposal and return
    return MarkovState(xi.parameters, xf_trajectory, compute_full_log_like(model, xi.parameters.value, xf_trajectory), prop_lk, prop_type)
end # end of std proposal function


## new standard proposal function
# NEED TO BENCHMARK AGAINST OLD ***
function standard_proposal(model::PrivateDiscuitModel, xi::MarkovState, xf_parameters::ParameterProposal)
    ## choose proposal type
    prop_type = rand(1:3)
    # trajectory proposal
    # - NEED TO MAKE THIS MORE EFFICIENT ****
    xf_trajectory = deepcopy(xi.trajectory)
    t0 = (model.t0_index == 0) ? 0.0 : xf_parameters.value[model.t0_index]
    if prop_type == 3
        ## move
        length(xi.trajectory.time) == 0 && (return MarkovState(xf_parameters, xi.trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, prop_type))
        # - IS THERE A MORE EFFICIENT WAY TO DO THIS? I.E. ROTATE using circshift or something?
        # choose event and define new one
        evt_i = rand(1:length(xi.trajectory.time))
        # evt_tm = (rand() * (model.obs_data.time[end] - t0)) + t0 #, xi.trajectory.event_type[evt_i])
        evt_tp = xi.trajectory.event_type[evt_i]
        # remove old one
        splice!(xf_trajectory.time, evt_i)
        splice!(xf_trajectory.event_type, evt_i)
        # add new one at random time
        add_event!(xf_trajectory, evt_tp, (rand() * (model.obs_data.time[end] - t0)) + t0)
        # compute ln g(x)
        prop_lk = 1.0
    else
        ## insert / delete
        if prop_type == 1
            ## choose type:
            tp = rand(1:size(model.m_transition, 1))
            ## insert at randomly chosen time
            add_event!(xf_trajectory, tp, (rand() * (model.obs_data.time[end] - t0)) + t0)
            ## compute ln g(x)
            prop_lk = (size(model.m_transition, 1) * (model.obs_data.time[end] - t0)) / length(xf_trajectory.time)
        else
            ## delete
            # println(" deleting... tp:", tp, " - ec: ", ec)
            length(xi.trajectory.time) == 0 && (return MarkovState(xi.parameters, xf_trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, prop_type))
            # choose event index (repeat if != tp)
            evt_i = rand(1:length(xi.trajectory.time))
            # remove
            splice!(xf_trajectory.time, evt_i)
            splice!(xf_trajectory.event_type, evt_i)
            # compute ln g(x)
            prop_lk = length(xi.trajectory.time) / ((model.obs_data.time[end] - t0) * size(model.m_transition, 1))
        end # end of insert/delete
    end
    ## evaluate full likelihood for trajectory proposal and return
    return MarkovState(xi.parameters, xf_trajectory, compute_full_log_like(model, xi.parameters.value, xf_trajectory), prop_lk, prop_type)
end # end of std proposal function

## model based proposal
include("discuit_mbp.jl")



## macros
# TBA? (proposal step)

## metropolis hastings algorithm (internal)
# - default proportion of parameter proposals (ppp): 0.3
# - NEED TO TEMPLATE FOR SINGLE EVENT TYPE MODELS ***********
function met_hastings_alg(model::PrivateDiscuitModel, steps::Int64, adapt_period::Int64, proposal_alg::Function, x0::MarkovState, prop_param::Bool, ppp::Float64) # full_like::Bool
    ## constants
    PARAMETER_PROPOSAL::Int64 = 4
    INITIAL_J::Float64 = 0.1
    # adaption interval
    a_h::Int64 = adapt_period / 10
    # initialise xi MOVE THIS OUT AND PASS ***********************
    xi = x0
    # covar matrix
    C_DEBUG && println( " initialising g for d = ", length(xi.parameters.value))
    covar = zeros(length(xi.parameters.value), length(xi.parameters.value))
    for i in eachindex(xi.parameters.value)
        covar[i,i] = 0.1 * xi.parameters.value[i] * xi.parameters.value[i]
    end
    g = Distributions.MvNormal(covar)
    sclr_j::Float64 = INITIAL_J
    # declare results
    mc = Array{Float64, 2}(undef, steps, length(xi.parameters.value))
    mcf = Array{Float64, 2}(undef, steps, length(xi.parameters.value))
    mc_log_like = Array{Float64,1}(undef, steps)
    # mc_prior = Array{Float64,1}(undef, steps)
    mc_accepted = falses(steps)
    # TO BE REMOVED? *******
    prop_type = Array{Int64,1}(undef, steps)
    ll_g = Array{Float64,1}(undef, steps)
    mh_p = Array{Float64,1}(undef, steps)
    mc_time = zeros(UInt64, steps)
    mc_mu = copy(xi.parameters.value)
    # is_tpd = exp(xi.parameters.prior + xi.log_like)
    # is_mu = is_tpd * xi.parameters.value

    # add first sample
    mc[1,:] .= xi.parameters.value
    mcf[1,:] .= xi.parameters.value
    mc_log_like[1] = xi.log_like
    # mc_prior[1] = xi.parameters.prior
    # mc_prior[1] = Distributions.pdf(xi.parameters.value)
    mc_accepted[1] = true
    # db
    prop_type[1] = 0
    ll_g[1] = 0
    mh_p[1] = 1
    st_time = time_ns()
    for i in 2:steps
        # make theta proposal
        adapt::Bool = true
        if prop_param
            # always make combined proposal (i.e. mbp)
            xf::MarkovState = proposal_alg(model, xi, get_mv_param(model, g, sclr_j, xi.parameters.value))
        else
            # choose parameter or trajectory proposal (i.e. standard)
            if rand() < ppp
                # parameter proposal
                prop = get_mv_param(model, g, sclr_j, xi.parameters.value)
                ll = Distributions.logpdf(model.prior, prop.value) == -Inf ? NULL_LOG_LIKE : compute_full_log_like(model, prop.value, xi.trajectory)
                xf = MarkovState(prop, xi.trajectory, ll, DF_PROP_LIKE, PARAMETER_PROPOSAL)
            else
                # trajectory proposal
                xf = proposal_alg(model, xi, xi.parameters)
                adapt = false
            end
        end
        # for DEBUG
        mcf[i,:] .= xf.parameters.value
        prop_type[i] = xf.prop_type
        ll_g[i] = xf.prop_like
        # PRIOR CHECK REDUNDANT? ******
        if (Distributions.logpdf(model.prior, xf.parameters.value) == -Inf || xf.log_like == NULL_LOG_LIKE)
            # reject automatically NEED TO THINK ABOUT THIS, BIT OF A MESS!
            mc_log_like[i] = NULL_LOG_LIKE
            mh_p[i] = -1.0
        else
            # accept or reject
            mc_log_like[i] = xf.log_like
            mh_prob::Float64 = exp(Distributions.logpdf(model.prior, xf.parameters.value) + xf.log_like - Distributions.logpdf(model.prior, xi.parameters.value) - xi.log_like)
            mh_p[i] = xf.prop_like * mh_prob
            if mh_prob > 1.0
                mc_accepted[i] = true
            else
                mh_prob > rand() && (mc_accepted[i] = true)
            end
        end
        ## acceptance handling
        if mc_accepted[i]
            mc[i,:] .= mcf[i,:]
            xi = xf
        else
            mc[i,:] .= mc[i - 1,:]
        end
        mc_time[i] = time_ns() - st_time
        if i < adapt_period
            ## ADAPTATION PERIOD
            # adjust theta jump scalar
            adapt && (sclr_j *= (mc_accepted[i] ? 1.002 : 0.999))
            # end of adaption period
            if i % a_h == 0
                # recalc covar matrix and update g
                covar = Distributions.cov(mc[1:i,:])
                if sum(covar) == 0
                    println("warning: low acceptance rate detected in adaptation period")
                else
                    # t0 stuffs
                    if model.t0_index > 0
                        covar[model.t0_index, model.t0_index] == 0.0 && (covar[model.t0_index, model.t0_index] = 1.0)
                    end
                    # update proposal dist
                    g = Distributions.MvNormal(covar)
                end
            end
        else
            ## ADAPTED
            # mu
            mc_mu .+=  mc[i,:]
            # is_mu .+= xf.parameters.value * exp(xf.parameters.prior + xf.log_like)
            # is_tpd += exp(xf.parameters.prior + xf.log_like)
        end
    end # end of Markov chain for loop
    # compute means and return results
    mc_mu ./= (steps - adapt_period)
    # is_mu ./= is_tpd
    pan = prop_param ? "MBP" : "Standard"
    ## MAKE GEWKE TEST OPTIONAL? ****************
    gw = run_geweke_test(mc, adapt_period)
    return MCMCResults(mc, mc_accepted, mc_mu, Distributions.cov(mc[(adapt_period + 1):steps,:]), pan, length(model.obs_data.time), adapt_period, gw, mcf, mc_log_like, xi, prop_type, ll_g, mh_p, mc_time)
end

## adaptive mh mcmc
# public wrappers
# - vanilla
function run_single_chain_analysis(model::DiscuitModel, obs_data::Observations, initial_parameters::Array{Float64, 1}, steps::Int64, adapt_period::Int64, mbp::Bool, ppp::Float64)
    # (model::DiscuitModel, obs_data::Observations, initial_parameters::Array{Float64, 1} = rand(model.prior), steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
    pm = get_private_model(model, obs_data)
    println("running single-chain ", mbp ? "MBP" : "std", "-MCMC analysis (model: ", model.model_name, ")")
    C_DEBUG && println(" theta init: ", initial_parameters)
    x0 = gillespie_sim_x0(pm, initial_parameters, !mbp)
    output = met_hastings_alg(pm, steps, adapt_period, mbp ? model_based_proposal : standard_proposal, x0, mbp, ppp)
    println(" finished (sample μ = ", output.mean, ")")
    return output
end
# - custom MH MCMC - called by multi chain function below
function run_custom_single_chain_analysis(model::DiscuitModel, obs_data::Observations, proposal_function::Function, x0::MarkovState, steps::Int64, adapt_period::Int64, prop_param::Bool, ppp::Float64)
# (model::DiscuitModel, obs_data::Observations, proposal_function::Function, x0::MarkovState, steps::Int64 = 50000, adapt_period::Int64 = 10000, prop_param::Bool = false, ppp::Float64 = 0.3)
    # ADD TIME / MSGS HERE *********************
    pm = get_private_model(model, obs_data)
    println("running custom-proposal MCMC analysis (model: ", model.model_name, ")")
    output =  met_hastings_alg(pm, steps, adapt_period, proposal_function, x0, prop_param, ppp)
    println(" finished (sample μ = ", output.mean, ").")
    return output
end

## convergence diagnostics
# autocorrelation R
"""
    compute_autocorrelation(mcmc, lags = 200)

**Parameters**
- `mcmc`    -- `MCMCResults` variable.
- `lags`    -- the number of lags to compute. Default: 200.

Compute autocorrelation R for a single Markov chain. Autocorrelation can be used to help determine how well the algorithm mixed by using `compute_autocorrelation(rs.mcmc)`. The autocorrelation function for a single Markov chain is implemented in Discuit using the standard formula:

\$R_l  = \\frac{\\textrm{E} [(X_i - \\bar{X})(X_{i+l} - \\bar{X})]}{\\sigma^2}\$

for any given lag `l` up to `lags` (default: 200).
"""
function compute_autocorrelation(mcmc::MCMCResults, lags::Int64 = 200)
    lag = zeros(Int64, lags + 1)
    output = zeros(lags + 1, length(mcmc.mean))
    # tmp = zeros(lags, length(mcmc.mean))
    # for each lag interval
    for l in 0:lags
        lag[l + 1] = l * AC_LAG_INT
        # for each parameter
        for j in 1:length(mcmc.mean)
            # compute SWAP AND VECTORISE
            for i in (mcmc.adapt_period + 1):(size(mcmc.samples, 1) - lag[l + 1])
                output[l + 1,j] += (mcmc.samples[i,j] - mcmc.mean[j]) * (mcmc.samples[i + lag[l + 1], j] - mcmc.mean[j])
                # tmp[l,j] += mcmc.samples[i,j]
            end
            output[l + 1,j] /= size(mcmc.samples, 1) - mcmc.adapt_period - lag[l + 1]
            output[l + 1,j] /= mcmc.covar[j,j]
            # tmp[l,j] /= size(mcmc.samples, 1) - mcmc.adapt_period - (l*AC_LAG_INT)
        end
    end
    # println("mu = ", mcmc.mean)
    return AutocorrelationResults(lag, output)
end
# autocorrelation R'
"""
    compute_autocorrelation(mcmc, lags = 200)

**Parameters**
- `mcmc`    -- an array of `MCMCResults` variables.
- `lags`    -- the number of lags to compute. Default: 200.

Compute autocorrelation R' for a two or more Markov chains. The formula for multiple chains is given by:

\$R^{\\prime}_l = \\frac{\\textrm{E} [ (X_i - \\bar{X}_b) ( X_{i + l} - \\bar{X}_b ) ]}{\\sigma^2_b}\$

\$\\sigma^2_b = \\textrm{E} [(X_i - \\bar{X}_b)^2]\$

for any given lag `l` up to `lags` (default: 200).
"""
function compute_autocorrelation(mcmc::Array{MCMCResults, 1}, lags::Int64 = 200)
    mce = Array{Float64, 2}(undef, length(mcmc), length(mcmc[1].mean))
    for mc in eachindex(mcmc)
        mce[mc,:] .= mcmc[mc].mean
    end
    # for each parameter
    mu = Array{Float64, 1}(undef, length(mcmc[1].mean))
    wcv = zeros(length(mcmc), length(mu))
    for j in eachindex(mu)
        # compute mean and whole chain var
        # - VECTORISE THIS? ********
        mu[j] = Statistics.mean(mce[:,j])
        for mc in eachindex(mcmc)
            for i in (mcmc[mc].adapt_period + 1):(size(mcmc[mc].samples, 1))
                wcv[j] += (mcmc[mc].samples[i,j] - mu[j]) * (mcmc[mc].samples[i,j] - mu[j])
            end
        end
        wcv[j] /= length(mcmc) * (size(mcmc[1].samples, 1) - mcmc[1].adapt_period)
    end
    # for each lag interval
    lag = zeros(Int64, lags + 1)
    output = zeros(lags + 1, length(mu))
    for l in 0:lags
        lag[l + 1] = l * AC_LAG_INT
        for j in eachindex(mu)
            for mc in eachindex(mcmc)
                for i in (mcmc[mc].adapt_period + 1):(size(mcmc[mc].samples, 1) - lag[l + 1])
                    output[l + 1,j] += (mcmc[mc].samples[i,j] - mu[j]) * (mcmc[mc].samples[i + lag[l + 1], j] - mu[j])
                end
            end
            output[l + 1,j] /= length(mcmc) * (size(mcmc[1].samples, 1) - mcmc[1].adapt_period)
            output[l + 1,j] /= wcv[j]
        end
    end
    return AutocorrelationResults(lag, output)
end

## Geweke test
function run_geweke_test(mc::Array{Float64,2}, adapt_period::Int64)
    NUM_A = 10
    ## compute interval
    n = size(mc, 1) - adapt_period
    h::Int64 = n / (NUM_A * 2)
    ## calculate b mean/variance
    nd2::Int64 = n / 2
    b_mu = Array{Float64, 1}(undef, size(mc, 2))
    b_var = Array{Float64, 1}(undef, size(mc, 2))
    for j in 1:size(mc, 2)
        b_mu[j] = Statistics.mean(mc[nd2:n,j])
        b_var[j] = Statistics.var(mc[nd2:n,j])
    end
    ## compute test statistics
    lbl = Array{Int64, 1}(undef, NUM_A)
    output = Array{Float64, 2}(undef, NUM_A, size(mc, 2))
    for i in 1:NUM_A
        frow = adapt_period + (i * h)
        lbl[i] = frow - adapt_period
        for j in 1:size(mc, 2)
            mu = Statistics.mean(mc[(frow - h):frow,j]) # collapse this ***************
            output[i,j] = (mu - b_mu[j]) / sqrt(Statistics.var(mc[(frow - h):frow,j]) + b_var[j])
        end
    end
    # return results as tuple
    return (lbl, output)
end

## enable multithreading
# NEED TO ADD DIFFERENT OS SUPPORT? ***
# function EnableMultithreading(number_of_threads::Int)
#     # linux, osx
#     export JULIA_NUM_THREADS=number_of_threads
#     # windows, c NEED TO ADD
# end

#### public wrappers ####

"""
    run_mcmc_analysis(model, obs_data, initial_parameters, steps = 50000, adapt_period = 10000, mbp = true, ppp = 0.3)

**Parameters**
- `model`               -- `DiscuitModel` (see [Discuit.jl models]@ref).
- `obs_data`            -- `Observations` data.
- `n_chains`            -- number of Markov chains (optional, default: 3.)
- `initial_parameters`  -- 2d array of initial model parameters. Each column vector correspondes to a single model parameter.
- `steps`               -- number of iterations.
- `adapt_period`        -- number of discarded samples.
- `mbp`                 -- model based proposals (MBP). Set `mbp = false` for standard proposals.
- `ppp`                 -- the proportion of parameter (vs. trajectory) proposals. Default: 30%. NB. not required for MBP.

Run an `n_chains`-MCMC analysis. The `initial_parameters` are sampled from the prior distribution unless otherwise specified by the user.

A Gelman-Rubin convergence diagnostic is automatically carried out for n_chains > 1 and included in the [multi-chain] results.

Otherwise the results of a single-chain analysis are returned, which include the Geweke test statistics computed for that analysis.
"""
function run_mcmc_analysis(model::DiscuitModel, obs_data::Observations; n_chains::Int64 = 3, initial_parameters = rand(model.prior, n_chains), steps::Int64 = 50000, adapt_period::Int64 = Int64(floor(steps / 5)), mbp::Bool = true, ppp::Float64 = 0.3)
    if size(initial_parameters, 2) == 1
        return run_single_chain_analysis(model, obs_data, initial_parameters[:,1], steps, adapt_period, mbp, ppp)
    else
        p_model = get_private_model(model, obs_data)
        ## initialise Markov chains
        println("running ", size(initial_parameters, 2) ,"-chain ", mbp ? "MBP" : "std", "-MCMC analysis (model: ", model.model_name, ")")
        C_DEBUG && println(" theta init: ", initial_parameters)
        mcmc = Array{MCMCResults,1}(undef, size(initial_parameters, 2))
        ## use @threads to multithread loop
        # Threads.@threads for i in eachindex(mcmc)
        for i in eachindex(mcmc)
            x0 = gillespie_sim_x0(p_model, initial_parameters[:,i], !mbp)
            mcmc[i] = met_hastings_alg(p_model, steps, adapt_period, mbp ? model_based_proposal : standard_proposal, x0, mbp, ppp)
            println(" chain ", i, " complete on thread ", Threads.threadid())
        end
        # REJOIN THREADS HERE ***
        ## ADD results check
        ## process results and return
        output = gelman_diagnostic(mcmc, size(initial_parameters, 1), steps - adapt_period)
        println(" finished (sample μ = ", round.(output.mu; sigdigits = C_PR_SIGDIG), ").")
        return output
    end
end

# - for custom MCMC
"""
    run_custom_mcmc_analysis(m_model, obs_data, proposal_function, x0, initial_parameters, steps = 50000, adapt_period = 10000, mbp = true, ppp = 0.3)

**Parameters**
- `model`               -- `DiscuitModel` (see [Discuit.jl models]@ref).
- `obs_data`            -- `Observations` data.
- `proposal_function`   -- `Function` for proposing changes to the trajectory. Must have the signature: `custom_proposal(model::PrivateDiscuitModel, xi::MarkovState, xf_parameters::ParameterProposal) = MarkovState(...)`
- `x0`                  -- vector of `MarkovState`s representing the initial samples.
- `initial_parameters`  -- matrix of initial model parameters. Each column vector correspondes to a single model parameter.
- `steps`               -- number of iterations.
- `adapt_period`        -- number of discarded samples.
- `prop_param`          -- simulaneously propose changes to parameters. Default: `false`.
- `ppp`                 -- the proportion of parameter (vs. trajectory) proposals. Default: 30%. NB. not required for MBP.

Run n (equal to the number of rows in `initial_parameters`) custom MCMC analyses and perform a Gelman-Rubin convergence diagnostic on the results. NEED TO OVERLOAD AND EXPAND.
"""
function run_custom_mcmc_analysis(model::DiscuitModel, obs_data::Observations, proposal_function::Function, x0::Array{MarkovState,1}; steps::Int64 = 50000, adapt_period::Int64 = Int64(floor(steps / 5)), prop_param::Bool = false, ppp::Float64 = 0.3)
# run_custom_single_chain_analysis
# (model::DiscuitModel, obs_data::Observations, proposal_function::Function, x0::MarkovState, steps::Int64 = 50000, adapt_period::Int64 = 10000, prop_param::Bool = false, ppp::Float64 = 0.3)
    if length(x0) == 1
        return run_custom_single_chain_analysis(model, obs_data, proposal_function, x0[1], steps, adapt_period, prop_param, ppp)
    else
        TEMP = 3
        model = get_private_model(model, obs_data)
        println("running ", length(x0), "-chain custom-proposal MCMC analysis (model: ", model.model_name, ")")
        ## initialise Markov chains
        mcmc = Array{MCMCResults,1}(undef, length(x0))
        ## CHANGE XO.PARAMS
        for i in eachindex(mcmc)
            mcmc[i] = met_hastings_alg(model, steps, adapt_period, proposal_function, x0[i], prop_param, ppp)
            println(" chain ", i, " complete.")
        end
        ## process results and return
        println(" - done. Running Gelman diagnostic...")
        output = gelman_diagnostic(mcmc, length(x0[1].parameters.value), steps - adapt_period)
        println(" finished (sample μ = ", output.mu, ").")
        return output
    end
end

## gelman diagnostic (internal)
function gelman_diagnostic(mcmc::Array{MCMCResults,1}, theta_size::Int64, num_iter::Int64)
    ## compute W; B; V
    # collect means and variances
    mce = Array{Float64, 2}(undef, length(mcmc), theta_size)
    mcv = Array{Float64, 2}(undef, length(mcmc), theta_size)
    # is_mu = zeros(theta_size)
    for i in eachindex(mcmc)
        mce[i,:] .= mcmc[i].mean
        # is_mu .+= mcmc[i].is_mu
        for j in 1:theta_size
            mcv[i,j] = mcmc[i].covar[j,j]
        end
    end
    # is_mu ./= length(mcmc)
    # compute W, B
    b = Array{Float64, 1}(undef, theta_size)
    w = Array{Float64, 1}(undef, theta_size)
    mu = Array{Float64, 1}(undef, theta_size)
    co = Array{Float64, 1}(undef, theta_size)
    v = Array{Float64, 1}(undef, theta_size)
    for j in 1:theta_size
        b[j] = num_iter * Distributions.cov(mce[:,j])
        w[j] = Statistics.mean(mcv[:,j])
        # mean of means and var of vars (used later)
        mu[j] = Statistics.mean(mce[:,j])
        co[j] = Distributions.cov(mcv[:,j])
        # compute pooled variance
        # - BENCHMARK THESE..?
        v[j] = w[j] * ((num_iter - 1) / num_iter) + b[j] * ((theta_size + 1) / (theta_size * num_iter))
        # v[j] = ((num_iter - 1) * (w[j] / num_iter)) + ((1 + (1 / size(theta_init, 1))) * (b[j] / num_iter))
        # v[j] = ((w[j] / num_iter) * (num_iter - 1)) + ((b[j] / num_iter) * (1 + (1 / size(theta_init, 1))))
    end
    #
    vv_w = Array{Float64, 1}(undef, theta_size)   # var of vars (theta_ex, i.e. W)
    vv_b = Array{Float64, 1}(undef, theta_size)   # var of vars (B)
    mce2 = copy(mce)
    mce2 .*= mce                                           # ev(theta)^2
    cv_wb = Array{Float64, 1}(undef, theta_size)   # wb covar
    for j in 1:theta_size
        vv_w[j] = co[j] / length(mcmc)
        vv_b[j] = (2 * b[j] * b[j]) / (length(mcmc) - 1)
        cv_wb[j] = (num_iter / length(mcmc)) * (Distributions.cov(mcv[:,j], mce2[:,j]) - (2 * mu[j] * Distributions.cov(mcv[:,j], mce[:,j])))
    end
    # compute d; d_adj (var.V)
    d = Array{Float64, 1}(undef, theta_size)
    dd = Array{Float64, 1}(undef, theta_size)
    atmp = num_iter - 1
    btmp = 1 + (1 / length(mcmc))
    for j in 1:theta_size
        tmp = ((vv_w[j] * atmp * atmp) + (vv_b[j] * btmp * btmp) + (cv_wb[j] * 2 * atmp * btmp)) / (num_iter * num_iter)
        d[j] = (2 * v[j] * v[j]) / tmp
        dd[j] = (d[j] + 3) / (d[j] + 1)
    end
    # compute scale reduction estimate
    sre = Array{Float64, 1}(undef, theta_size)
    sre_ll = Array{Float64, 1}(undef, theta_size)
    sre_ul = Array{Float64, 1}(undef, theta_size)
    for j in 1:theta_size
        rr = (1 + (1 / length(mcmc))) * (1 / num_iter)  * (b[j] / w[j])
        sre[j] = sqrt(dd[j] * (((num_iter - 1) / num_iter) + rr))
        # F dist(nu1, nu2)
        fdst = Distributions.FDist(length(mcmc) - 1, 2 * w[j] * w[j] / vv_w[j])
        sre_ll[j] = sqrt(dd[j] * (((num_iter - 1) / num_iter) + Distributions.quantile(fdst, 0.025) * rr))
        sre_ul[j] = sqrt(dd[j] * (((num_iter - 1) / num_iter) + Distributions.quantile(fdst, 0.975) * rr))
    end
    # return results
    return MultiMCMCResults(mu, sqrt.(w), sre, sre_ll, sre_ul, mcmc)
end

include("./discuit_utils.jl")

## MBP IBIS - need to fold this in at some point ***
include("discuit_ibis.jl")

end # end of module
