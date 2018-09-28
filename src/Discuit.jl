"""
# module Discuit

- Julia version:
- Author: martin
- Date: 2018-08-22

# Examples

```jldoctest
julia>
```
"""
module Discuit

## depends
using Random
using Distributions
using CSV

## exports
# public structs
export DiscuitModel, ObsData, SimResults
# core functionality
export set_random_seed, gillespie_sim, run_met_hastings_mcmc, compute_autocorrelation, run_gelman_diagnostic
# utilities (e.g. saving results to file0
export print_trajectory, print_observations, print_mcmc_results, print_autocorrelation, print_gelman_results, read_obs_data_from_file
# custom functionality (in development)
export MarkovState, ParameterProposal, Event, PrivateDiscuitModel, run_custom_mcmc, run_custom_mcmc_gelman_diagnostic
# for unit testing:
export PrivateDiscuitModel, get_private_model, model_based_proposal, standard_proposal, gillespie_sim_x0, run_geweke_test, compute_full_log_like


include("./discuit_structs.jl")

## constants
# max trajectory constant (make this optional?)
const MAX_TRAJ = 196000
const DF_PROP_LIKE = 1.0
const NULL_LOG_LIKE = -Inf

## model generating functions
# TO BE ADDED *******

## FOR DEBUG
"""
    set_random_seed(seed)

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
function iterate_sim!(model::PrivateDiscuitModel, trajectory::Array{Event, 1}, population::Array{Int32,1}, parameters::Array{Float64,1}, time::Float64, tmax::Float64)
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
        push!(trajectory, Event(time, et))
        # update population
        population .+= model.m_transition[et,:] # event_types[et].v_transition
        # print("\n t: ", time, ". R: ", cum_rates[end])
    end
end
## NEED TO TEMPLATE FOR OPTIONAL EVENTS
# gillespie_sim
# gillespie_sim_x0
# MBP
# Std proposal

# run sim and return trajectory
# - ADD option for >1 sim?
"""
    gillespie_sim(model, parameters, tmax = 100.0, num_obs = 5)

Run a DGA simulation on `model`. Returns a SimResults containing the trajectory and observations data.
"""
function gillespie_sim(model::DiscuitModel, parameters::Array{Float64,1}, tmax::Float64 = 100.0, num_obs::Int64 = 5)
    # initialise some things
    p_model = get_private_model(model, ObsData([0.0], [0 0]) ) # NOTE TO JAMIE: can this be made nullable? or is it better to use delegates?
    obs_times = collect(tmax / num_obs : tmax / num_obs : tmax)
    obs_vals = Array{Int32, 2}(undef, length(obs_times), length(p_model.initial_condition))
    # time = 0.0
    population = copy(p_model.initial_condition)
    trajectory = Event[]
    # run
    t_prev = model.t0_index == 0 ? 0.0 : parameters[model.t0_index]
    for i in eachindex(obs_times)
        iterate_sim!(p_model, trajectory, population, parameters, t_prev, obs_times[i])
        obs_vals[i,:] .= p_model.obs_function(population)
        t_prev = obs_times[i]
    end
    # return trajectory
    return SimResults(trajectory, ObsData(obs_times, obs_vals))
end
# sim to initialise Markov chain
function gillespie_sim_x0(model::PrivateDiscuitModel, parameters::Array{Float64,1}, full_like::Bool)
    generate = true
    while generate
        # initialise some things
        population = copy(model.initial_condition)
        trajectory = Event[]
        output = 0.0    # obs model
        # run
        t_prev = model.t0_index == 0 ? 0.0 : parameters[model.t0_index]
        # print("\n t0: ", t_prev)
        for i in eachindex(model.obs_data.time)
            # print("\n period: ", i, ". time: ", t_prev)
            iterate_sim!(model, trajectory, population, parameters, t_prev, model.obs_data.time[i])
            # evaluate log likelihood
            output += model.obs_model(model.obs_data.val[i,:], population)
            output == NULL_LOG_LIKE && break
            t_prev = model.obs_data.time[i]
        end
        ## REPLACE WITH CONST OR PARAMETER? *****
        if output != NULL_LOG_LIKE && length(trajectory) > 5
            return MarkovState(ParameterProposal(parameters, model.prior(parameters)), trajectory, full_like ? compute_full_log_like(model, parameters, trajectory) : output, DF_PROP_LIKE, 0)
        end
    end
end
## mv parameter proposal
function get_mv_param(model::PrivateDiscuitModel, g::MvNormal, sclr::Float64, theta_i::Array{Float64, 1})
    output = rand(g)
    output .*= sclr
    output .+= theta_i
    return ParameterProposal(output, model.prior(output))
end
## standard trajectory proposal
# likelihood function
function compute_full_log_like(model::PrivateDiscuitModel, parameters::Array{Float64,1}, trajectory::Array{Event, 1})
    population = copy(model.initial_condition)
    # log_like = 0.0
    ll_traj = 0.0
    ll_obs = 0.0
    t = model.t0_index == 0 ? 0.0 : parameters[model.t0_index]
    evt_i = 1
    # workspace
    lambda = Array{Float64, 1}(undef, size(model.m_transition, 1))
    # for each observation period:
    for obs_i in eachindex(model.obs_data.time)
        while evt_i <= length(trajectory)
            trajectory[evt_i].time > model.obs_data.time[obs_i] && break
            model.rate_function(lambda, parameters, population)
            ## deal with event
            # event log likelihood
            ll_traj += log(lambda[trajectory[evt_i].event_type]) - (sum(lambda) * (trajectory[evt_i].time - t))
            # evt_i < 4 && print("\n i: ", evt_i, ". ll = ln ", lambda[trajectory[evt_i].event_type], " - ", sum(lambda), "*", trajectory[evt_i].time - t ,"=", ll_traj)
            # update population and handle -ve (NB. template? check not req'd for MBP)
            population .+= model.m_transition[trajectory[evt_i].event_type,:]
            any(x->x<0, population) && return -Inf
            t = trajectory[evt_i].time
            evt_i += 1
        end
        ## handle observation
        model.rate_function(lambda, parameters, population)
        # log likelihood
        # model.pop_index > 0 && (ll_traj += log(lambda[model.pop_index]))
        ll_traj -= sum(lambda) * (model.obs_data.time[obs_i] - t)
        # obs model
        ll_obs += model.obs_model(model.obs_data.val[obs_i,:], population)
        ll_traj == NULL_LOG_LIKE && return log(0.0)
        t = model.obs_data.time[obs_i]
    end
    # model.fc_index > 0 && (population[model.fc_index] == 0 || (log_like = NULL_LOG_LIKE))
    # print("\n ll: ", ll_traj, " ::: ", ll_obs)
    return ll_traj + ll_obs
end
# event type count
function get_event_type_count(trajectory::Array{Event,1}, et::Int64)
    output = 0
    for i in eachindex(trajectory)
        trajectory[i].event_type == et && (output += 1)
    end
    return output
end
# proposal function
function standard_proposal(model::PrivateDiscuitModel, xi::MarkovState, xf_parameters::ParameterProposal)
    ## choose proposal type
    prop_type = rand(1:3)
    # trajectory proposal
    # - NEED TO MAKE THIS MORE EFFICIENT ****
    xf_trajectory = copy(xi.trajectory)
    t0 = (model.t0_index == 0) ? 0.0 : theta_f.value[model.t0_index]
    if prop_type == 3
        ## move
        length(xi.trajectory) == 0 && (return MarkovState(theta_f, xi.trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, prop_type))
        # - IS THERE A MORE EFFICIENT WAY TO DO THIS? I.E. ROTATE using circshift or something?
        # choose event and define new one
        evt_i = rand(1:length(xi.trajectory))
        evt = Event((rand() * (model.obs_data.time[end] - t0)) + t0, xi.trajectory[evt_i].event_type)
        # remove old one
        splice!(xf_trajectory, evt_i)
        # add new one
        if evt.time > xf_trajectory[end].time
            push!(xf_trajectory, evt)
        else
            for i in eachindex(xf_trajectory)
                if xf_trajectory[i].time > evt.time
                    insert!(xf_trajectory, i, evt)
                    break
                end
            end
        end
        # compute ln g(x)
        prop_lk = 1.0
    else
        ## insert / delete
        # choose type and count
        tp = rand(1:size(model.m_transition, 1))
        ec = get_event_type_count(xf_trajectory, tp)
        if prop_type == 1
            ## insert
            # choose time
            tm = (rand() * (model.obs_data.time[end] - t0)) + t0
            # insert at new index
            if (length(xf_trajectory) == 0 || tm > xf_trajectory[end].time)
                push!(xf_trajectory, Event(tm, tp))
            else
                for i in eachindex(xf_trajectory)
                    if xf_trajectory[i].time > tm
                        insert!(xf_trajectory, i, Event(tm, tp))
                        break
                    end
                end
            end
            # compute ln g(x)
            prop_lk = (model.obs_data.time[end] - t0) / (ec + 1)
        else
            ## delete
            # print("\n deleting... tp:", tp, " - ec: ", ec)
            ec == 0 && (return MarkovState(xi.parameters, xf_trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, prop_type))
            # choose event index (repeat if != tp)
            evt_i = rand(1:length(xi.trajectory))
            while xi.trajectory[evt_i].event_type != tp
                evt_i = rand(1:length(xi.trajectory))
            end
            # remove
            splice!(xf_trajectory, evt_i)
            # compute ln g(x)
            prop_lk = ec / (model.obs_data.time[end] - t0)
        end # end of insert/delete
    end
    ## evaluate full likelihood for trajectory proposal and return
    return MarkovState(xi.parameters, xf_trajectory, compute_full_log_like(model, xi.parameters.value, xf_trajectory), prop_lk, prop_type)
end # end of std proposal function
## model based proposal
# single mbp iteration
function iterate_mbp(model::PrivateDiscuitModel, obs_i::Int64, evt_i::Int64, time::Float64, xi::MarkovState, pop_i::Array{Int32, 1}, xf_trajectory::Array{Event, 1}, theta_f::Array{Float64, 1}, pop_f::Array{Int32, 1})
    # workspace
    lambda_i = Array{Float64, 1}(undef, size(model.m_transition, 1))
    lambda_f = Array{Float64, 1}(undef, size(model.m_transition, 1))
    lambda_d = Array{Float64, 1}(undef, size(model.m_transition, 1))
    # iterate until next observation
    while true
        if evt_i > length(xi.trajectory)
            tmax = model.obs_data.time[obs_i]
        else
            tmax = min(model.obs_data.time[obs_i], xi.trajectory[evt_i].time)
        end
        model.rate_function(lambda_i, xi.parameters.value, pop_i)
        while true
            # calculate rate delta
            model.rate_function(lambda_f, theta_f, pop_f)
            for i in 1:size(model.m_transition, 1)
                lambda_d[i] = max(lambda_f[i] - lambda_i[i], 0.0)
            end
            cumsum!(lambda_d, lambda_d)
            # 0 rate test
            lambda_d[end] == 0.0 && break
            time -= log(rand()) / lambda_d[end]
            # break if event time exceeded
            time > tmax && break
            # else choose event type (init as final event)
            et = choose_event(lambda_d)
            # add event to trajectory
            push!(xf_trajectory, Event(time[1], et))
            length(xf_trajectory) > MAX_TRAJ && (return evt_i)
            # update population
            pop_f .+= model.m_transition[et,:]
            # t_max check
        end
        ## handle event
        # checks
        evt_i > length(xi.trajectory) && break
        xi.trajectory[evt_i].time > model.obs_data.time[obs_i] && break
        # else
        et = xi.trajectory[evt_i].event_type
        time = xi.trajectory[evt_i].time
        # prob_keep::Float64 = get_p_keep(et
        prob_keep = lambda_f[et] / lambda_i[et]
        keep = false
        if prob_keep > 1.0
            keep = true
        else
            prob_keep > rand() && (keep = true)
        end
        if keep
            # add to trajectory and update population
            push!(xf_trajectory, Event(time, et))
            pop_f .+= model.m_transition[et,:]
        end
        # update for i regardless
        pop_i .+= model.m_transition[et,:]
        # - iterate event counter
        evt_i += 1
    end
    return evt_i
end
# mbp function
function model_based_proposal(model::PrivateDiscuitModel, xi::MarkovState, xf_parameters::ParameterProposal)
    MBP_PROP_TYPE = 10
    # make theta proposal (RENAME TO T F)
    # prop = get_mv_param(model, g, sclr, theta_i.value)
    if xf_parameters.prior == 0.0
        # no need to evaluate
        return MarkovState(xf_parameters, xi.trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, MBP_PROP_TYPE)
    else
        # initialise
        xf_trajectory = Event[]
        pop_i = copy(model.initial_condition)
        pop_f = copy(model.initial_condition)
        # - time / event counter
        time = model.t0_index == 0 ? 0.0 : prop.value[model.t0_index]
        evt_i = [1]
        # - log likelihood
        output = 0.0
        # for each observation period
        for obs_i in eachindex(model.obs_data.time)
            # iterate until next observation
            evt_i[1] = iterate_mbp(model, obs_i, evt_i[1], time, xi, pop_i, xf_trajectory, xf_parameters.value, pop_f)
            # check trajectory length
            length(xf_trajectory) > MAX_TRAJ && (return MarkovState(xf_parameters, xi.trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, MBP_PROP_TYPE))
            # else handle observation
            time = model.obs_data.time[obs_i]
            output += model.obs_model(model.obs_data.val[obs_i,:], pop_f)
        end
        return MarkovState(xf_parameters, xf_trajectory, output, DF_PROP_LIKE, MBP_PROP_TYPE)
    end
end

## adaptive mh mcmc
# public wrappers
# - vanilla
"""
    run_met_hastings_mcmc(model, obs_data, initial_parameters, steps = 50000, adapt_period = 10000, mbp = true, ppp = 0.3)

Run an MCMC analysis based on `model` and `obs_data` of type ObsData. The number of samples obtained will be `steps` - `adapt_period`. Set `mbp` = true for model based proposals (false for standard). The proportion of parameter proposals (`ppp`) can be changed from the default value of 30%.
"""
function run_met_hastings_mcmc(model::DiscuitModel, obs_data::ObsData, initial_parameters::Array{Float64, 1}, steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
    # ADD TIME / MSGS HERE *********************
    pm = get_private_model(model, obs_data)
    print("\nrunning MCMC.")
    x0 = gillespie_sim_x0(pm, initial_parameters, !mbp)
    # print(string("\nx0 length: ", length(x0.trajectory)))
    return met_hastings_alg(pm, steps, adapt_period, mbp ? model_based_proposal : standard_proposal, x0, mbp, ppp)
end
# - custom MH MCMC
"""
    run_custom_mcmc(model, obs_data, proposal_function, x0, steps = 50000, adapt_period = 10000, prop_param = false, ppp = 0.3)

Run a custom MCMC analysis. Similar to `run_met_hastings_mcmc` except that the`proposal_function` (of type Function) and initial state `x0` (of type MarkovState) are user defined.
"""
function run_custom_mcmc(model::DiscuitModel, obs_data::ObsData, proposal_function::Function, x0::MarkovState, steps::Int64 = 50000, adapt_period::Int64 = 10000, prop_param::Bool = false, ppp::Float64 = 0.3)
    # ADD TIME / MSGS HERE *********************
    pm = get_private_model(model, obs_data)
    print("\nrunning custom MCMC.")
    return met_hastings_alg(pm, steps, adapt_period, proposal_function, x0, prop_param, ppp)
end

# metropolis hastings algorithm
# - default proportion of parameter proposals (ppp): 0.3
function met_hastings_alg(model::PrivateDiscuitModel, steps::Int64, adapt_period::Int64, proposal_alg::Function, x0::MarkovState, prop_param::Bool, ppp::Float64) # full_like::Bool
    ## constants
    PARAMETER_PROPOSAL::Int64 = 4
    INITIAL_J::Float64 = 0.1
    # adaption interval
    a_h::Int64 = adapt_period / 10
    # initialise xi MOVE THIS OUT AND PASS ***********************
    xi = x0
    # covar matrix
    covar = zeros(length(xi.parameters.value), length(xi.parameters.value))
    for i in eachindex(xi.parameters.value)
        covar[i,i] = 0.1 * xi.parameters.value[i] * xi.parameters.value[i]
    end
    g = MvNormal(covar)
    sclr_j::Float64 = INITIAL_J
    # declare results
    mc = Array{Float64, 2}(undef, steps, length(xi.parameters.value))
    mcf = Array{Float64, 2}(undef, steps, length(xi.parameters.value))
    mc_log_like = Array{Float64,1}(undef, steps)
    mc_prior = Array{Float64,1}(undef, steps)
    mc_accepted = falses(steps)
    # TO BE REMOVED? *******
    prop_type = Array{Int64,1}(undef, steps)
    ll_g = Array{Float64,1}(undef, steps)
    mh_p = Array{Float64,1}(undef, steps)
    # get some samples
    # - NEED TO TIDY THIS UP ***
    mc[1,:] .= xi.parameters.value
    mcf[1,:] .= xi.parameters.value
    mc_log_like[1] = xi.log_like
    mc_prior[1] = xi.parameters.prior
    mc_accepted[1] = true
    # db
    prop_type[1] = 0
    ll_g[1] = 0
    mh_p[1] = 1
    for i in 2:steps
        ## JP, CAN THIS BE DONE BETTER WITH TEMPLATE?
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
                ll = prop.prior == 0.0 ? NULL_LOG_LIKE : compute_full_log_like(model, prop.value, xi.trajectory)
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
        # mc_prior[i] = xf.parameters.prior
        # PRIOR CHECK REDUNDANT? ******
        if (xf.parameters.prior == 0.0 || xf.log_like == NULL_LOG_LIKE)
            # reject automatically NEED TO THINK ABOUT THIS, BIT OF A MESS!
            mc_log_like[i] = NULL_LOG_LIKE
            mh_p[i] = -1.0
        else
            mc_log_like[i] = xf.log_like
            # accept or reject
            mh_prob::Float64 = xf.prop_like * (xf.parameters.prior / xi.parameters.prior) * exp(xf.log_like - xi.log_like)
            mh_p[i] = mh_prob
            if mh_prob > 1.0
                mc_accepted[i] = true
            else
                mh_prob > rand() && (mc_accepted[i] = true)
            end
        end
        #
        if mc_accepted[i]
            mc[i,:] .= mcf[i,:]
            xi = xf
        else
            mc[i,:] .= mc[i - 1,:]
        end
        ## ADAPTATION PERIOD
        if i < adapt_period
            # adjust theta jump scalar
            adapt && (sclr_j *= (mc_accepted[i] ? 1.002 : 0.999))
            # end of adaption period
            if i % a_h == 0
                # recalc covar matrix and update g
                covar = cov(mc[1:i,:])
                if sum(covar) == 0
                    print("\nwarning: low acceptance rate detected in adaptation period")
                else
                    # t0 stuffs
                    if model.t0_index > 0
                        covar[model.t0_index, model.t0_index] == 0.0 && (covar[model.t0_index, model.t0_index] = 1.0)
                    end
                    # print("\n covar: ")
                    # print(covar)
                    g = MvNormal(covar)
                end
            end
        end
    end # end of Markov chain for loop
    # compute mean/var and return results
    mc_bar = Array{Float64, 1}(undef, length(xi.parameters.value))
    for i in eachindex(mc_bar)
        mc_bar[i] = mean(mc[(adapt_period + 1):steps,i])
    end
    pan = prop_param ? "MBP" : "Standard"
    ## MAKE GEWKE TEST OPTIONAL ****************
    gw = run_geweke_test(mc, adapt_period)
    return MCMCResults(mc, mc_accepted, mc_bar, cov(mc[(adapt_period + 1):steps,:]), pan, length(model.obs_data.time), adapt_period, gw, mcf, mc_log_like, xi, prop_type, ll_g, mh_p)
end
## convergence diagnostics
const LAG_INT = 10
# autocorrelation R
function compute_autocorrelation(mcmc::MCMCResults, lags::Int64 = 200)
    output = zeros(lags + 1, length(mcmc.mean))
    # tmp = zeros(lags, length(mcmc.mean))
    # for each lag interval
    for l in 0:lags
        # for each parameter
        for j in 1:length(mcmc.mean)
            # compute SWAP AND VECTORISE
            for i in (mcmc.adapt_period + 1):(size(mcmc.samples, 1) - (l*LAG_INT))
                output[l + 1,j] += (mcmc.samples[i,j] - mcmc.mean[j]) * (mcmc.samples[i + (l*LAG_INT), j] - mcmc.mean[j])
                # tmp[l,j] += mcmc.samples[i,j]
            end
            output[l + 1,j] /= size(mcmc.samples, 1) - mcmc.adapt_period - (l*LAG_INT)
            output[l + 1,j] /= mcmc.covar[j,j]
            # tmp[l,j] /= size(mcmc.samples, 1) - mcmc.adapt_period - (l*LAG_INT)
        end
    end
    # print(tmp[1,:])
    print("\nmu = ", mcmc.mean)

    return output
end
# autocorrelation R'
# - ADD OPTION FOR STANDARD FORM..? *
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
        mu[j] = mean(mce[:,j])
        for mc in eachindex(mcmc)
            for i in (mcmc[mc].adapt_period + 1):(size(mcmc[mc].samples, 1))
                wcv[j] += (mcmc[mc].samples[i,j] - mu[j]) * (mcmc[mc].samples[i,j] - mu[j])
            end
        end
        wcv[j] /= length(mcmc) * (size(mcmc[1].samples, 1) - mcmc[1].adapt_period)
    end
    # for each lag interval
    output = zeros(lags + 1, length(mu))
    for l in 0:lags
        for j in eachindex(mu)
            for mc in eachindex(mcmc)
                for i in (mcmc[mc].adapt_period + 1):(size(mcmc[mc].samples, 1) - (l*LAG_INT))
                    output[l + 1,j] += (mcmc[mc].samples[i,j] - mu[j]) * (mcmc[mc].samples[i + (l*LAG_INT), j] - mu[j])
                end
            end
            output[l + 1,j] /= length(mcmc) * (size(mcmc[1].samples, 1) - mcmc[1].adapt_period)
            output[l + 1,j] /= wcv[j]
        end
    end
    return output
end

## Geweke test
# NB. would it be better to do this for the whole chain? **********
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
        b_mu[j] = mean(mc[nd2:n,j])
        b_var[j] = var(mc[nd2:n,j])
    end
    ## compute test statistics
    lbl = Array{Int64, 1}(undef, NUM_A)
    output = Array{Float64, 2}(undef, NUM_A, size(mc, 2))
    for i in 1:NUM_A
        frow = adapt_period + (i * h)
        lbl[i] = frow - adapt_period
        for j in 1:size(mc, 2)
            mu = mean(mc[(frow - h):frow,j]) # collapse this ***************
            output[i,j] = (mu - b_mu[j]) / sqrt(var(mc[(frow - h):frow,j]) + b_var[j])
        end
    end
    # return results as tuple
    return (lbl, output)
end

## Gelman-Rubin diagnostic
# NEED TO ADD OVERLOAD for those already run
# public wrappers
# - calls mcmc
"""
    run_gelman_diagnostic(m_model, obs_data, initial_parameters, steps = 50000, adapt_period = 10000, mbp = true, ppp = 0.3)

Run n (equal to the number of rows in `initial_parameters`)  MCMC analyses and perform a Gelman-Rubin convergence diagnostic on the results. NEED TO OVERLOAD AND EXPAND.
"""
function run_gelman_diagnostic(m_model::DiscuitModel, obs_data::ObsData, initial_parameters::Array{Float64, 2}, steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
    model = get_private_model(m_model, obs_data)
    print("\nrunning gelman diagnostic...")
    ## initialise Markov chains
    # NEED TO ADD MULTI THREADING OR ASYNC HERE ***************
    mcmc = Array{MCMCResults,1}(undef, size(initial_parameters, 1))
    for i in eachindex(mcmc)
        x0 = gillespie_sim_x0(model, initial_parameters[i,:], !mbp)
        # print(string("\nx0 length: ", length(x0.trajectory)))
        mcmc[i] = met_hastings_alg(model, steps, adapt_period, mbp ? model_based_proposal : standard_proposal, x0, mbp, ppp)
        print("\n chain ", i, " complete")
    end
    # REJOIN THREADS HERE **********************
    ## ADD results check
    ## process results and return
    return gelman_diagnostic(mcmc, size(initial_parameters, 2), steps - adapt_period)
end
# for custom MCMC
# - TO BE REMOVED? ************
function run_custom_mcmc_gelman_diagnostic(model::DiscuitModel, obs_data::ObsData, proposal_function::Function, x0::Array{MarkovState,1}, steps::Int64 = 50000, adapt_period::Int64 = 10000, prop_param::Bool = false, ppp::Float64 = 0.3)
    TEMP = 3
    model = get_private_model(model, obs_data)
    print("\nrunning gelman diagnostic...")
    ## initialise Markov chains
    mcmc = Array{MCMCResults,1}(undef, length(x0))
    ## CHANGE XO.PARAMS
    for i in eachindex(mcmc)
        mcmc[i] = met_hastings_alg(model, steps, adapt_period, proposal_function, x0[i], prop_param, ppp)
        print("\n chain ", i, " complete")
    end
    ## process results and return
    return gelman_diagnostic(mcmc, length(x0[1].parameters.value), steps - adapt_period)
end
# internal function
function gelman_diagnostic(mcmc::Array{MCMCResults,1}, theta_size::Int64, num_iter::Int64)
    ## compute W; B; V
    # collect means and variances
    mce = Array{Float64, 2}(undef, length(mcmc), theta_size)
    mcv = Array{Float64, 2}(undef, length(mcmc), theta_size)
    for i in eachindex(mcmc)
        mce[i,:] .= mcmc[i].mean
        for j in 1:theta_size
            mcv[i,j] = mcmc[i].covar[j,j]
        end
    end
    # compute W, B
    b = Array{Float64, 1}(undef, theta_size)
    w = Array{Float64, 1}(undef, theta_size)
    mu = Array{Float64, 1}(undef, theta_size)
    co = Array{Float64, 1}(undef, theta_size)
    v = Array{Float64, 1}(undef, theta_size)
    for j in 1:theta_size
        b[j] = num_iter * cov(mce[:,j])
        w[j] = mean(mcv[:,j])
        # mean of means and var of vars (used later)
        mu[j] = mean(mce[:,j])
        co[j] = cov(mcv[:,j])
        # compute pooled variance
        v = Array{Float64, 1}(undef, theta_size)
    end
    # compute pooled variance
    # v = Array{Float64, 1}(undef, size(theta_init, 2))
    # for j in 1:size(theta_init, 2)
    #     # BENCHMARK THESE..?
    #     v[j] = w[j] * ((num_iter - 1) / num_iter) + b[j] * ((size(theta_init, 1) + 1) / (size(theta_init, 1) * num_iter))
    #     # v[j] = ((num_iter - 1) * (w[j] / num_iter)) + ((1 + (1 / size(theta_init, 1))) * (b[j] / num_iter))
    #     # v[j] = ((w[j] / num_iter) * (num_iter - 1)) + ((b[j] / num_iter) * (1 + (1 / size(theta_init, 1))))
    # end
    #
    vv_w = Array{Float64, 1}(undef, theta_size)   # var of vars (theta_ex, i.e. W)
    vv_b = Array{Float64, 1}(undef, theta_size)   # var of vars (B)
    mce2 = copy(mce)
    mce2 .*= mce                                           # ev(theta)^2
    cv_wb = Array{Float64, 1}(undef, theta_size)   # wb covar
    for j in 1:theta_size
        vv_w[j] = co[j] / length(mcmc)
        vv_b[j] = (2 * b[j] * b[j]) / (length(mcmc) - 1)
        cv_wb[j] = (num_iter / length(mcmc)) * (cov(mcv[:,j], mce2[:,j]) - (2 * mu[j] * cov(mcv[:,j], mce[:,j])))
    end
    # compute d; d_adj (var.V)
    # - SHOULD BENCHMARK PRECOMP **********
    d = Array{Float64, 1}(undef, theta_size)
    dd = Array{Float64, 1}(undef, theta_size)
    # print("\n var.v: ")
    atmp = num_iter - 1
    btmp = 1 + (1 / length(mcmc))
    for j in 1:theta_size
        tmp = ((vv_w[j] * atmp * atmp) + (vv_b[j] * btmp * btmp) + (cv_wb[j] * 2 * atmp * btmp)) / (num_iter * num_iter)
        # print(" ", tmp)
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
        fdst = FDist(length(mcmc) - 1, 2 * w[j] * w[j] / vv_w[j])
        sre_ll[j] = sqrt(dd[j] * (((num_iter - 1) / num_iter) + quantile(fdst, 0.025) * rr))
        sre_ul[j] = sqrt(dd[j] * (((num_iter - 1) / num_iter) + quantile(fdst, 0.975) * rr))
    end
    # return results
    return GelmanResults(mu, sre, sre_ll, sre_ul, mcmc)
end
## print convergence test results
# print Geweke results (private)

# print autocorrelation
"""
    print_autocorrelation(autocorrelation::Array{Float64, 2}, fpath::String)

Save the results from a call to compute_autocorrelation to the file `fpath`, e.g. "./out/ac.csv".
"""
function print_autocorrelation(autocorrelation::Array{Float64, 2}, fpath::String)
    open(fpath, "w") do f
        # print headers
        write(f, "lag")
        for j in 1:size(autocorrelation, 2)
            write(f, ", x$j")
        end
        # print autocorr
        for i in 1:size(autocorrelation, 1)
            write(f, "\n$((i - 1) * LAG_INT)")
            for j in 1:size(autocorrelation, 2)
                write(f, ", $(autocorrelation[i,j])")
            end
        end
    end
end
# print Gelman results
"""
    print_gelman_results(results::GelmanResults, dpath::String)

Save the results from a call to run_gelman_diagnostic to the directory `dpath`, e.g. "./out/gelman/".
"""
function print_gelman_results(results::GelmanResults, dpath::String)
    # create directory
    # dpath = string("./out/", dname, "/")
    isdir(dpath) || mkpath(dpath)
    # print summary by theta row
    open(string(dpath, "gelman.csv"), "w") do f
        # print headers
        write(f, "theta, mu, sre, sre_ll, sre_ul")
        for p in eachindex(results.mu)
            write(f, "\n$p, $(results.mu[p]), $(results.sre[p]), $(results.sre_ll[p]), $(results.sre_ul[p])")
        end
    end # end of print summary
    # print chains
    for i in eachindex(results.mcmc)
        print_mcmc_results(results.mcmc[i], string(dpath, "mc", i, "/"))
    end
end
## print MCMC results to file
"""
    print_mcmc_results(mcmc::MCMCResults, dpath::String)

Save the results from a call to run_met_hastings_mcmc or run_custom_mcmc to the directory `dpath`, e.g. "./out/mcmc/".
"""
function print_mcmc_results(mcmc::MCMCResults, dpath::String)
    # NEED TO ADD / IF NOT THERE ALREADY *******
    # create directory
    isdir(dpath) || mkpath(dpath)
    # print metadata
    open(string(dpath, "metadata.csv"), "w") do f
        # print headers
        write(f, "proposal_alg, num_obs, adapt_period")
        # print md
        write(f, "\n$(mcmc.proposal_alg), $(mcmc.num_obs), $(mcmc.adapt_period)")
    end
    # print theta summary
    open(string(dpath, "theta.csv"), "w") do f
        # print headers
        write(f, "theta, mean, sd")
        for p in eachindex(mcmc.mean)
            sd = sqrt(mcmc.covar[p,p])
            write(f, "\n$p, $(mcmc.mean[p]), $sd")
        end
    end
    # print Geweke results
    open(string(dpath, "geweke.csv"), "w") do f
        # print headers
        write(f, "x")
        for p in eachindex(mcmc.mean)
            write(f, ", $p")
        end
        # print test statistics
        for i in eachindex(mcmc.geweke[1])
            write(f, "\n$(mcmc.geweke[1][i])")
            for p in eachindex(mcmc.mean)
                write(f, ", $(mcmc.geweke[2][i,p])")
            end
        end
    end
    # print samples
    open(string(dpath, "chain.csv"), "w") do f
        # print headers
        write(f, "iter, accepted")
        for p in 1:size(mcmc.samples, 2)
            write(f, ", x$p, xf$p")
        end
        write(f, ", xf_ll, prop_type, ll_g, mh_prob")
        #
        for i in 1:size(mcmc.samples, 1)
            write(f, "\n $i, $(mcmc.mc_accepted[i])")
            for p in 1:size(mcmc.samples, 2)
                # t = mc[i, p]
                write(f, ", $(mcmc.samples[i, p]), $(mcmc.mcf[i, p])")
            end
            write(f, ", $(mcmc.mc_log_like[i]), $(mcmc.prop_type[i]), $(mcmc.ll_g[i]), $(mcmc.mh_prob[i])")
        end
    end # end of print
end
## save trajectory to file
"""
    print_trajectory(model, sim_results, fpath)

Save an augmented trajectory from a variable of type `SimResults` to the file `fpath`, e.g. "./out/sim.csv".
"""
function print_trajectory(model::DiscuitModel, sim_results::SimResults, fpath::String)
    open(fpath, "w") do f
        population = copy(model.initial_condition)
        # print headers
        write(f, "time, event")
        for p in eachindex(population)
            c = model.model_name[p]
            write(f, ", $c")
        end
        # print event trajectory
        evt_i = 1
        for obs_i in eachindex(sim_results.observations.time)
            # handle events
            while evt_i <= length(sim_results.trajectory)
                sim_results.trajectory[evt_i].time > sim_results.observations.time[obs_i] && break
                tp = sim_results.trajectory[evt_i].event_type
                write(f, "\n $(sim.trajectory[evt_i].time), $tp")
                population .+= model.m_transition[tp,:]
                for p in 1:length(population)
                    write(f, ", $(population[p])")
                end
                evt_i += 1
            end
            # handle observation
            model.obs_function(population)
            write(f, "\n $(sim.observations.time[obs_i]), -1")
            for p in 1:length(population)
                write(f, ", $(population[p])")
            end
        end
    end # end of print
end # end of function
## save observations data to file
"""
    print_observations(obs_data, fpath)

Save a set of observations (e.g. from a `SimResults` variable) to the file `fpath`, e.g. "./out/obs.csv".
"""
function print_observations(obs_data::ObsData, fpath::String)
    open(fpath, "w") do f
        # print headers
        write(f, "time")
        for j in 1:size(obs_data.val, 2)
            write(f, ", val$j")
        end
        # print observations
        for i in eachindex(obs_data.time)
            write(f, "\n$(obs_data.time[i])")
            for j in 1:size(obs_data.val, 2)
                write(f, ", $(obs_data.val[i,j])")
            end
        end
    end # end of print
end
"""
    read_obs_data_from_file(fpath)

Read a set of observations from the location `fpath` and return the results as a variable of type `ObsData`.
"""
function read_obs_data_from_file(fpath::String)
    # NEED TO FIX WARNING ***
    df = CSV.read(fpath)
    return ObsData(df[1], df[:,2:size(df, 2)])
end

end # end of module