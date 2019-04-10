### model based proposal
# this code is based on the method described by Pooley et al. 2015

## single mbp iteration
function iterate_mbp(model::PrivateDiscuitModel, obs_i::Int64, evt_i::Int64, time::Float64, xi::MarkovState, pop_i::Array{Int64}, xf_trajectory::Trajectory, theta_f::Array{Float64, 1}, pop_f::Array{Int64})
    # workspace
    lambda_i = Array{Float64, 1}(undef, size(model.m_transition, 1))
    lambda_f = Array{Float64, 1}(undef, size(model.m_transition, 1))
    lambda_d = Array{Float64, 1}(undef, size(model.m_transition, 1))
    # iterate until next observation
    while true
        if evt_i > length(xi.trajectory.time)
            tmax = model.obs_data.time[obs_i]
        else
            tmax = min(model.obs_data.time[obs_i], xi.trajectory.time[evt_i])
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
            push!(xf_trajectory.time, time)
            push!(xf_trajectory.event_type,  et)
            length(xf_trajectory.time) > MAX_TRAJ && (return evt_i)
            # update population
            pop_f .+= model.m_transition[et,:]
        end
        ## handle event
        # checks
        evt_i > length(xi.trajectory.time) && break
        xi.trajectory.time[evt_i] > model.obs_data.time[obs_i] && break
        # else
        et = xi.trajectory.event_type[evt_i]
        time = xi.trajectory.time[evt_i]
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
            push!(xf_trajectory.time, time)
            push!(xf_trajectory.event_type, et)
            pop_f .+= model.m_transition[et,:]
        end
        # update for i regardless
        pop_i .+= model.m_transition[et,:]
        # - iterate event counter
        evt_i += 1
    end
    return evt_i
end

## initialise sequence for MBP
function initialise_sequence(model::PrivateDiscuitModel, xi::MarkovState, pop_i::Array{Int64}, xf_trajectory::Trajectory, theta_f::Array{Float64, 1}, pop_f::Array{Int64})
    evt_i::Int64 = 1
    if theta_f[model.t0_index] < xi.parameters.value[model.t0_index]
        # workspace
        lambda_f = Array{Float64, 1}(undef, size(model.m_transition, 1))
        # sim on 'full'
        time::Float64 = theta_f[model.t0_index]
        while true
            # calculate rate delta
            model.rate_function(lambda_f, theta_f, pop_f)
            cumsum!(lambda_f, lambda_f)
            # 0 rate test
            lambda_f[end] == 0.0 && break
            time -= log(rand()) / lambda_f[end]
            # break if prev t0 time exceeded
            time > xi.parameters.value[model.t0_index] && break
            # else choose event type (init as final event)
            et = choose_event(lambda_f)
            # add event to trajectory
            push!(xf_trajectory.time, time)
            push!(xf_trajectory.event_type,  et)
            length(xf_trajectory.time) > MAX_TRAJ && (return evt_i)
            # update population
            pop_f .+= model.m_transition[et,:]
        end
    else
        # 'delete'
        while true
            evt_i > length(xi.trajectory.time) && break
            xi.trajectory.time[evt_i] > theta_f[model.t0_index] && break
            # update for i
            pop_i .+= model.m_transition[xi.trajectory.event_type[evt_i],:]
            # - iterate event counter
            evt_i += 1
        end
    end
    return evt_i
end

## mbp function
function model_based_proposal(model::PrivateDiscuitModel, xi::MarkovState, xf_parameters::ParameterProposal)
    MBP_PROP_TYPE = 10
    # make theta proposal (RENAME TO T F)
    if xf_parameters.prior == 0.0
        # no need to evaluate
        return MarkovState(xf_parameters, xi.trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, MBP_PROP_TYPE)
    else
        # initialise
        xf_trajectory = Trajectory(Float64[], Int64[])
        pop_i = copy(model.initial_condition)
        pop_f = copy(model.initial_condition)
        # - time / event counter
        evt_i = [1]
        time = 0.0
        if model.t0_index > 0
            evt_i[1] = initialise_sequence(model, xi, pop_i, xf_trajectory, xf_parameters.value, pop_f)
            time = max(xf_parameters.value[model.t0_index], xi.parameters.value[model.t0_index])
        end
        # - log likelihood
        output = 0.0
        # for each observation period
        for obs_i in eachindex(model.obs_data.time)
            # iterate until next observation
            evt_i[1] = iterate_mbp(model, obs_i, evt_i[1], time, xi, pop_i, xf_trajectory, xf_parameters.value, pop_f)
            # check trajectory length
            length(xf_trajectory.time) > MAX_TRAJ && (return MarkovState(xf_parameters, xi.trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, MBP_PROP_TYPE))
            # else handle observation
            time = model.obs_data.time[obs_i]
            output += model.observation_model(model.obs_data.val[obs_i,:], pop_f)
        end
        return MarkovState(xf_parameters, xf_trajectory, output, DF_PROP_LIKE, MBP_PROP_TYPE)
    end
end
