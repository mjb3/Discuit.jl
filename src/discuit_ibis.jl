##################### RESOURCES #####################
import LinearAlgebra
import Statistics
import PrettyTables

##################### STRUCTS #####################

## event
struct Event
    time::Float64
    event_type::Int64
end

## a single realisation of the model
struct Particle
    theta::Array{Float64,1}
    initial_condition::Array{Int64}
    final_condition::Array{Int64}
    trajectory::Array{Event,1}
    log_like::Array{Float64,1}  # prior, g(x)
end

## IBIS sample
struct ImportanceSample
    mu::Array{Float64,1}
    cv::Array{Float64,2}
    theta::Array{Float64,2}
    weight::Array{Float64,1}
    run_time::UInt64
    bme::Array{Float64,1}
end

##################### RESAMPLE #####################
## systematic (samples single seed u(0,1/N])
# Carpenter (1999)
function rs_systematic(w::Array{Float64,1})
    cw = cumsum(w)
    output = Array{Int64,1}(undef, length(w))
    u = Array{Float64,1}(undef, length(w))
    u[1] = rand() / length(w) # sample ~ U(0,1/N]
    for i in 2:length(cw)
        u[i] = u[1] + ((i - 1) / length(w))
    end
    u .*= cw[end]
    # set output = new index
    j = 1
    for i in eachindex(output)
        while u[i] > cw[j]
            j += 1
        end
        output[i] = j
    end
    return output
end

##################### DGA #####################

## gillespie sim iteration
# - NEED TO REDO WITH MACROS ***************
function iterate_particle!(p::Particle, model::PrivateDiscuitModel, time::Float64, yi::Int64) #tmax::Float64
    cum_rates = Array{Float64, 1}(undef, size(model.m_transition,1))
    while true
        model.rate_function(cum_rates, p.theta, p.final_condition)
        cumsum!(cum_rates, cum_rates)
        cum_rates[end] == 0.0 && break          # 0 rate test
        time -= log(rand()) / cum_rates[end]
        time > model.obs_data.time[yi] && break                  # break if max time exceeded
        et = choose_event(cum_rates)            # else choose event type (init as final event)
        p.final_condition .+= model.m_transition[et,:]  # update population CHANGED ***
        push!(p.trajectory, Event(time, et))    # add event to sequence
        if length(p.trajectory) > MAX_TRAJ      # HACK
            p.log_like[2] = -Inf
            return p.log_like[2]
        end
    end
    output = model.observation_model(model.obs_data.val[yi,:], p.final_condition)
    # output = model.observation_model(y, p.final_condition, p.theta)
    p.log_like[2] += output
    return output
end

##################### CMN #####################

## Gaussian mv parameter proposal
function get_mv_param(propd::Distributions.MvNormal, sclr::Float64, theta_i::Array{Float64, 1})
    output = rand(propd)
    output .*= sclr
    output .+= theta_i
    return output
end

##################### MBP #####################

## iterator
function iterate_mbp!(xf, pop_i::Array{Int64}, model::PrivateDiscuitModel, obs_i::Int64, evt_i::Int64, time::Float64, xi)
    ## workspace
    # - zeros(model.n_events)
    lambda_f = Array{Float64,1}(undef, size(model.m_transition,1))
    lambda_i = Array{Float64,1}(undef, size(model.m_transition,1))
    lambda_d = Array{Float64,1}(undef, size(model.m_transition,1))
    ## iterate until next observation
    while true
        tmax = evt_i > length(xi.trajectory) ? model.obs_data.time[obs_i] : min(model.obs_data.time[obs_i], xi.trajectory[evt_i].time)
        model.rate_function(lambda_i, xi.theta, pop_i)
        while true
            model.rate_function(lambda_f, xf.theta, xf.final_condition)
            for i in eachindex(lambda_d)            # compute rate delta
                lambda_d[i] = max(lambda_f[i] - lambda_i[i], 0.0)
            end
            cumsum!(lambda_d, lambda_d)
            lambda_d[end] == 0.0 && break           # 0 rate test
            time -= log(rand()) / lambda_d[end]
            time > tmax && break                    # break if event time exceeded
            et = choose_event(lambda_d)             # else choose event type (init as final event)
            xf.final_condition .+= model.m_transition[et,:]
            push!(xf.trajectory, Event(time, et))   # add event to trajectory
        end
        ## handle event
        evt_i > length(xi.trajectory) && break
        xi.trajectory[evt_i].time > model.obs_data.time[obs_i] && break
        et::Int64 = xi.trajectory[evt_i].event_type
        time = xi.trajectory[evt_i].time
        prob_keep = lambda_f[et] / lambda_i[et]
        if (prob_keep > 1.0 || prob_keep > rand())
            push!(xf.trajectory, Event(time, et))
            xf.final_condition .+= model.m_transition[et,:]
        end
        pop_i .+= model.m_transition[et,:]
        evt_i += 1
    end
    return evt_i
end

## initialise sequence for MBP
function initialise_trajectory!(xf, pop_i::Array{Int64}, model::PrivateDiscuitModel, xi)
    # xf_trajectory::Array{Event,1}
    evt_i::Int64 = 1
    if model.t0_index == 0
        return evt_i
    else
        if xf.theta[model.t0_index] < xi.theta[model.t0_index]
            ## 'sim'
            lambda_f = zeros(size(model.m_transition,1))                    # workspace
            t = xf.theta[model.t0_index]
            while true
                # compute rates:
                model.rate_function(lambda_f, xf.theta, xf.final_condition)
                cumsum!(lambda_f, lambda_f)
                lambda_f[end] == 0.0 && break                   # 0 rate test
                t -= log(rand()) / lambda_f[end]
                t > xi.theta[model.t0_index] && break           # break if prev t0 time exceeded
                et = choose_event(lambda_f)                     # else choose event type (init as final event)
                push!(xf.trajectory, Event(t, et))              # add event to trajectory
                xf.final_condition .+= model.m_transition[et,:] # update population
                length(xf.trajectory) > MAX_TRAJ && (return evt_i)
            end
        else
            ## 'delete'
            while true
                evt_i > length(xi.trajectory) && break
                xi.trajectory[evt_i].time > xf.theta[model.t0_index] && break
                pop_i .+= model.m_transition[xi.trajectory[evt_i].event_type,:] # update population
                evt_i += 1                                                      # iterate event counter
            end
        end
        return evt_i
    end
end

## up to time y_i
function partial_model_based_proposal(model::PrivateDiscuitModel, theta_f::Array{Float64,1}, xi::Particle, ymax::Int64)
    # xf = Particle(theta_f, copy(xi.initial_condition), copy(xi.initial_condition), Event[], [log(model.prior_density(theta_f)), 0.0, 0.0])
    xf = Particle(theta_f, copy(xi.initial_condition), copy(xi.initial_condition), Event[], [Distributions.logpdf(model.prior, theta_f), 0.0, 0.0])
    ## evaluate prior density
    if xf.log_like[1] == -Inf
        xf.log_like[2] = xf.log_like[1]
        return xf
    else
        # debug && println("init ll: ", xf.log_like[2])
        ## make mbp
        pop_i = copy(xi.initial_condition) # GET RID *
        evt_i = initialise_trajectory!(xf, pop_i, model, xi)
        time::Float64 = model.t0_index == 0 ? 0.0 : max(xf.theta[model.t0_index], xi.theta[model.t0_index])
        for obs_i in 1:ymax
            ## iterate MBP
            evt_i = iterate_mbp!(xf, pop_i, model, obs_i, evt_i, time, xi)
            if length(xf.trajectory) > MAX_TRAJ # - HACK: MAX trajectory check
                xf.log_like[2] = -Inf
                return xf
            end
            time = model.obs_data.time[obs_i]   ## handle observation
            xf.log_like[3] = model.observation_model(model.obs_data.val[obs_i,:], xf.final_condition)
            xf.log_like[2] += xf.log_like[3]
        end
        return xf
    end
end

##################### UTILS #####################

## compute sd (internal)
function compute_sigma(cv::Array{Float64,2})
    sd = zeros(size(cv,1))
    for i in eachindex(sd)
        sd[i] = sqrt(cv[i,i])
    end
    return sd
end

## importance sample:
function tabulate_results(results::ImportanceSample, proposals = false)
    ## samples
    println("IBIS results:")
    h = ["θ", "μ", "σ", "BME"]
    d = Matrix(undef, length(results.mu), 4)
    sd = compute_sigma(results.cv)
    d[:,1] .= 1:length(results.mu)
    d[:,2] .= round.(results.mu; sigdigits = C_PR_SIGDIG)
    d[:,3] .= round.(sd; sigdigits = C_PR_SIGDIG)
    d[:,4] .= 0
    d[1:2,4] = round.(results.bme; sigdigits = C_PR_SIGDIG)
    PrettyTables.pretty_table(d, h)
end

## print theta summary (internal use)
function print_sample_summary(results, dpath::String)
    # print theta summary
    open(string(dpath, "summary.csv"), "w") do f
        # print headers
        write(f, "theta,mu,sigma")
        # print data
        for p in eachindex(results.mu)
            # c = model.model_name[p]
            write(f, "\n$p,$(results.mu[p]),$(sqrt(results.cv[p,p]))")
        end
    end
end

## print importance sample (just the weighed sample and summary)
function print_imp_sample(results::ImportanceSample, dpath::String)
    # check dir
    isdir(dpath) || mkpath(dpath)
    # print importance samples
    open(string(dpath, "theta.csv"), "w") do f
        # print headers
        write(f, "1")
        for i in 2:length(results.mu)
            write(f, ",$i")
        end
        # print data
        for i in 1:size(results.theta, 1)
            write(f, "\n$(results.theta[i,1])")
            for p in 2:length(results.mu)
                write(f, ",$(results.theta[i,p])")
            end
        end
    end
    # print weights
    open(string(dpath, "weight.csv"), "w") do f
        # print headers
        write(f, "i,w")
        for i in eachindex(results.weight)
            write(f, "\n$i,$(results.weight[i])")
        end
    end
    # print theta summary
    print_sample_summary(results, string(dpath, "is_"))
end

"""
    print_results

**Parameters**
- `samples`     -- a data structure of type `MCMCSample` or `ImportanceSample`.
- `dpath`       -- the directory where the results will be saved.

Print the results of an inference analysis to file.
"""
## print importance sample results
function print_results(results::ImportanceSample, dpath::String)
    # check dir
    isdir(dpath) || mkpath(dpath)
    # print metadata
    open(string(dpath, "metadata.csv"), "w") do f
        write(f, "st,n,run_time,bme\nis,$(length(results.mu)),$(results.run_time),$(results.bme[1])")
    end
    ##
    print_imp_sample(results, dpath)
end

##################### MAIN #####################

## compute effective sample size
function compute_ess(w::Array{Float64,1})
    return sum(w)^2 / sum(w.^2)
end

## compute is mu var
function compute_is_mu_covar!(mu::Array{Float64,1}, cv::Array{Float64,2}, theta::Array{Float64,2}, w::Array{Float64,1})
    for i in eachindex(mu)
        mu[i] = sum(w .* theta[:,i]) / sum(w)                       # compute mu
        cv[i,i] = sum(w .* ((theta[:,i] .- mu[i]).^2)) / sum(w)     # variance
        for j in 1:(i-1)                                            # covar
            cv[i,j] = cv[j,i] = sum(w .* (theta[:,i] .- mu[i]) .* (theta[:,j] .- mu[j])) / sum(w)
        end
    end
end

##
function get_prop_density(cv::Array{Float64,2}, old)
    ## update proposal density
    tmp = LinearAlgebra.Hermitian(cv)
    if LinearAlgebra.isposdef(tmp)
        return Distributions.MvNormal(Matrix(tmp))
    else
        println("warning: particle degeneracy problem")
        println(" covariance matrix: ", cv)
        return old
    end
end

## MBP IBIS algorithm
C_DF_ESS_CRIT = 0.5
"""
    run_mbp_ibis_analysis(model, obs_data, initial_parameters, ess_rs_crit = 0.5; n_props = 3, ind_prop = false, alpha = 1.002)

**Parameters**
- `model`               -- `DiscuitModel` (see [Discuit.jl models]@ref).
- `obs_data`            -- `Observations` data.
- `np`                  -- number of particles (default = 2000.)
- `ess_rs_crit`         -- resampling criteria (default = 0.5.)
- `n_props`             -- MBP mutations per step (default = 3.)
- `ind_prop`            -- true for independent theta proposals (default = false.)
- `alpha`               -- user-defined, increase for lower acceptance rate targeted (default = 1.002.)

Run an MBP IBIS analysis based on `model` and `obs_data` of type `Observations`, with manually specified initial theta.
"""
function run_mbp_ibis_analysis(mdl::DiscuitModel, obs_data::Observations; np = 2000, ess_rs_crit = C_DF_ESS_CRIT, n_props = 3, ind_prop = false, alpha = 1.002)
    # , theta_init::Array{Float64, 2}
    theta_init = transpose(rand(mdl.prior, np))
    ## initialise
    model = get_private_model(mdl, obs_data)
    outer_p = np # TIDY THIS UP *********************
    println("running MBP IBIS analysis for n = ", outer_p, " (model: ", mdl.model_name, ")")
    start_time = time_ns()
    ess_crit = ess_rs_crit * outer_p
    # fn_rs = rsp_systematic

    ## initialise particles
    ptcls = Array{Particle,1}(undef, outer_p)
    theta = copy(theta_init) # GET RID? *
    for p in eachindex(ptcls)
        ic = model.initial_condition
        # ptcls[p] = Particle(theta[p,:], ic, copy(ic), Event[], [log(model.prior_density(theta[p,:])), 0.0])
        ptcls[p] = Particle(theta[p,:], ic, copy(ic), Event[], [Distributions.logpdf(model.prior, theta[p,:]), 0.0, 0.0])
    end
    ## resampling workspace
    ptcls2 = deepcopy(ptcls)

    ## proposal distribution and scalar
    propd = Distributions.MvNormal(Matrix{Float64}(LinearAlgebra.I, size(theta_init,2), size(theta_init,2)))
    tj = 0.2

    ## initialise and run:
    w = ones(outer_p)
    mu = zeros(size(theta, 2))
    cv = zeros(size(theta, 2), size(theta, 2))
    k_log = zeros(Int64, 2)
    ## for each observation
    bme = zeros(2)
    gx = zeros(outer_p)
    mtd_gx = zeros(outer_p)
    t = zeros(outer_p)
    model.t0_index > 0 && (t .= theta[:, model.t0_index])
    for obs_i in eachindex(model.obs_data.time)
        ## for each 'outer' particle
        for p in eachindex(ptcls)
            ## compute incremental weights (i.e. run pf)
            gx[p] = exp(iterate_particle!(ptcls[p], model, t[p], obs_i))
        end
        ## COMPUTE L and update weights
        lml = log(sum(w .* gx) / sum(w))
        bme[1] += lml
        # bme[1] += log(sum(w .* gx) / sum(w))
        w .*= gx
        ##
        compute_is_mu_covar!(mu, cv, theta, w)
        ## resample and mutate if criteria satisfied:
        essv = compute_ess(w)
        # t = model.obs_data[obs_i].time
        if (essv < ess_crit)
            ## resample and swap
            propd = get_prop_density(cv, propd)
            nidx = rs_systematic(w)
            # println(obs_i, " - ", length(nidx), " - ", length(mtd_gx))
            mtd_gx .= gx[nidx]
            for p in eachindex(ptcls)
                ptcls2[p] = deepcopy(ptcls[nidx[p]])
            end
            mlr = Statistics.mean(gx[nidx]) * exp(lml)
            ptcls, ptcls2 = ptcls2, ptcls
            # mutate:
            k_log[1] += outer_p * n_props
            for p in eachindex(ptcls)
                for mki in 1:n_props
                    ## propose new theta - independent OR conditional on current sample (rec'd)
                    theta_f = ind_prop ? get_mv_param(propd, 1.0, mu) : get_mv_param(propd, tj, ptcls[p].theta)
                    # - RETRIEVE MOST RECENT MARGINAL
                    xf = partial_model_based_proposal(model, theta_f, ptcls[p], obs_i)
                    ## HACK: ADD BACK?
                    if exp(sum(xf.log_like) - sum(ptcls[p].log_like)) > rand()
                    # if exp(xf.log_like[2] - ptcls[p].log_like[2]) > rand()
                        mtd_gx[p] = exp(xf.log_like[3])
                        ptcls[p] = xf
                        k_log[2] += 1
                        tj *= alpha
                    else
                        tj *= 0.999
                    end
                end
                theta[p,:] .= ptcls[p].theta
            end
            ## RB ML update
            bme[2] += log(mlr / Statistics.mean(mtd_gx))
            w .= 1  # reset w = 1
        else
            ## standard ML update
            bme[2] += log(sum(w .* gx) / sum(w))
        end # END OF sample/mutate
        # obs_min = obs_i + 1
        t .= model.obs_data.time[obs_i]
    end # END OF OBS LOOP
    compute_is_mu_covar!(mu, cv, theta, w)
    println(" - IS mu = ", mu)
    # bme[1] += log(sum(w) / outer_p)
    # bme[2] += log(Statistics.mean(w))
    bme .*= -2

    ## return results
    println(" - finished. AR = ", round(100.0 * k_log[2] / k_log[1]; sigdigits = 3), "%")
    return ImportanceSample(mu, cv, theta, w, time_ns() - start_time, bme)#, rejs
end

## draw initial samples
function draw_from_limit(n::Int64, theta_limits::Array{Float64, 1})
    output = zeros(n, length(theta_limits))
    for i in eachindex(theta_limits)
        output[:,i] .= rand(n) * theta_limits[i]
    end
    return output
end
