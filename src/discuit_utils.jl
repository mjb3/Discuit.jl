## print autocorrelation
"""
    print_autocorrelation(autocorrelation, fpath)

**Parameters**
- `autocorrelation` -- the results of a call to `compute_autocorrelation`.
- `fpath`           -- the file path of the destination file.

Save the results from a call to `compute_autocorrelation` to the file `fpath`, e.g. "./out/ac.csv".
"""
function print_autocorrelation(acr::AutocorrelationResults, fpath::String)
    open(fpath, "w") do f
        # print headers
        write(f, "lag")
        for j in 1:size(acr.autocorrelation, 2)
            write(f, ",x$j")
        end
        # print autocorr
        for i in 1:size(acr.autocorrelation, 1)
            # write(f, "\n$((i - 1) * AC_LAG_INT)")
            write(f, "\n$(acr.lag[i])")
            for j in 1:size(acr.autocorrelation, 2)
                write(f, ",$(acr.autocorrelation[i,j])")
            end
        end
    end
end

## print Gelman results
"""
    print_results(results::MultiMCMCResults, dpath::String)

**Parameters**
- `results` -- `MultiMCMCResults` variable.
- `dpath`   -- the path of the directory where the results will be saved.

Save the results from a call to `run_multi_chain_analysis` to the directory `dpath`, e.g. "./out/gelman/".
"""
function print_results(results::MultiMCMCResults, dpath::String)
    # create directory
    # dpath = string("./out/", dname, "/")
    isdir(dpath) || mkpath(dpath)
    # run metadata
    open(string(dpath, "metadata.csv"), "w") do f
        write(f, "chains\n$(length(results.mcmc))")
    end
    # print summary by theta row
    open(string(dpath, "parameters.csv"), "w") do f
        # print headers
        write(f, "parameter,mu,sdw,sre,sre_ll,sre_ul")
        for p in eachindex(results.mu)
            write(f, "\n$p,$(results.mu[p]),$(results.sdw[p]),$(results.sre[p]),$(results.sre_ll[p]),$(results.sre_ul[p])")
        end
    end # end of print summary
    # print chains
    for i in eachindex(results.mcmc)
        print_results(results.mcmc[i], string(dpath, "mc", i, "/"))
    end
end

## print MCMC results to file
"""
    print_results(mcmc, dpath)

**Parameters**
- `results` -- `MCMCResults` variable.
- `dpath`   -- the path of the directory where the results will be saved.

Save the results from a call to `run_single_chain_analysis` or `run_custom_single_chain_analysis` to the directory `dpath`, e.g. "./out/mcmc/".
"""
function print_results(mcmc::MCMCResults, dpath::String)
    # NEED TO ADD / IF NOT THERE ALREADY *******
    # create directory
    isdir(dpath) || mkpath(dpath)
    # print metadata
    open(string(dpath, "metadata.csv"), "w") do f
        # print headers
        write(f, "proposal_alg,num_obs,adapt_period")
        # print md
        write(f, "\n$(mcmc.proposal_alg),$(mcmc.num_obs),$(mcmc.adapt_period)")
    end
    # print parameter summary
    open(string(dpath, "parameters.csv"), "w") do f
        # print headers
        write(f, "parameter,mean,sd")
        for p in eachindex(mcmc.mean)
            sd = sqrt(mcmc.covar[p,p])
            write(f, "\n$p,$(mcmc.mean[p]),$sd")
        end
    end
    # print Geweke results
    open(string(dpath, "geweke.csv"), "w") do f
        # print headers
        write(f, "lag")
        for p in eachindex(mcmc.mean)
            write(f, ",$p")
        end
        # print test statistics
        for i in eachindex(mcmc.geweke[1])
            write(f, "\n$(mcmc.geweke[1][i])")
            for p in eachindex(mcmc.mean)
                write(f, ",$(mcmc.geweke[2][i,p])")
            end
        end
    end
    # print samples
    open(string(dpath, "samples.csv"), "w") do f
        # print headers
        write(f, "iter,accepted")
        for p in 1:size(mcmc.samples, 2)
            write(f, ",x$p,xf$p")
        end
        write(f, ",xf_ll,prop_type,ll_g,mh_prob,sys_time")
        #
        for i in 1:size(mcmc.samples, 1)
            write(f, "\n $i,$(mcmc.mc_accepted[i])")
            for p in 1:size(mcmc.samples, 2)
                # t = mc[i, p]
                write(f, ",$(mcmc.samples[i, p]),$(mcmc.mcf[i, p])")
            end
            write(f, ",$(mcmc.mc_log_like[i]),$(mcmc.prop_type[i]),$(mcmc.ll_g[i]),$(mcmc.mh_prob[i]),$(mcmc.mc_time[i])")
        end
    end # end of print
end

## save trajectory to file
"""
    print_trajectory(model, sim_results, fpath)

**Parameters**
- `model`       -- `DiscuitModel` (see [Discuit.jl models]@ref).
- `sim_results` -- `SimResults` variable.
- `fpath`       -- the destination file path.

Save an augmented trajectory from a variable of type `SimResults` (i.e. from a call to `gillespie_sim`, see [Simulation](@ref)) to the file `fpath`, e.g. "./out/sim.csv".
"""
function print_trajectory(model::DiscuitModel, sim_results::SimResults, fpath::String)
    ndims(model.m_transition) > 2 && throw("sorry, can't handle population dim >1")
    # check dir
    # isdir(dpath) || mkpath(dpath)
    open(fpath, "w") do f
        population = copy(model.initial_condition)
        # print headers
        write(f, "time,event")
        for p in eachindex(population)
            c = model.model_name[p]
            write(f, ",$c")
        end
        # print event trajectory
        evt_i = 1
        for obs_i in eachindex(sim_results.observations.time)
            # handle events
            while evt_i <= length(sim_results.trajectory.time)
                sim_results.trajectory.time[evt_i] > sim_results.observations.time[obs_i] && break
                tp = sim_results.trajectory.event_type[evt_i]
                write(f, "\n $(sim_results.trajectory.time[evt_i]),$tp")
                population .+= model.m_transition[tp,:]
                for p in 1:length(population)
                    write(f, ",$(population[p])")
                end
                evt_i += 1
            end
            # handle observation
            model.observation_function(population)
            write(f, "\n $(sim_results.observations.time[obs_i]), -1")
            for p in 1:length(population)
                write(f, ",$(population[p])")
            end
        end
    end # end of print
end # end of function

## save observations data to file
"""
    print_observations(obs_data, fpath)

**Parameters**
- `obs_data`    -- `Observations` data.
- `fpath`       -- the destination file path.

Save a set of observations (e.g. from a `SimResults` obtained by a call to `gillespie_sim` to the file `fpath`, e.g. "./out/obs.csv".
"""
function print_observations(obs_data::Observations, fpath::String)
    ndims(obs_data.val) > 2 && throw("sorry, can't handle population dim >2")
    open(fpath, "w") do f
        # print headers
        write(f, "time")
        for j in 1:size(obs_data.val, 2)
            write(f, ",val$j")
        end
        # print observations
        for i in eachindex(obs_data.time)
            write(f, "\n$(obs_data.time[i])")
            for j in 1:size(obs_data.val, 2)
                write(f, ",$(obs_data.val[i,j])")
            end
        end
    end # end of print
end

## get observations data from object or file location
"""
    get_observations(source)

**Parameters**
- `source`      -- `Array`, `DataFrame` or filepath (i.e. `String`) containing the data (with times in the first column).

Create and return a variable of type `Observations` based on a two dimensional array, `DataFrame` or file location.
"""
function get_observations(array::Array{Float64, 2})
    y = Array{Int64, 2}(undef, size(array, 1), size(array, 2) - 1)
    y .= array[:,2:size(array, 2)]
    return Observations(array[:,1], y)
end
function get_observations(df::DataFrames.DataFrame)
    return Observations(df[1], df[2:size(df, 2)])
end
function get_observations(fpath::String)
    df = CSV.read(fpath)
    return get_observations(df)
end

### tabulate stuff
import PrettyTables
const C_PR_SIGDIG = 3

## MCMC proposal summary
function tabulate_proposals(results::MCMCResults)
    println("Proposal summary:")
    h = ["Adapted", "Proposed", "Accepted", "Rate"]
    n_iter = length(results.mc_accepted)
    adp_prd = results.adapt_period
    ## summarise:
    pr = [adp_prd, n_iter - adp_prd]
    ac = zeros(Int64, 2)
    ac[1] = sum(results.mc_accepted[1:adp_prd])
    ac[2] = sum(results.mc_accepted[(adp_prd + 1):n_iter])
    d = Matrix(undef, 3, 4)
    d[1, :] .= [ "false", pr[1], ac[1], round(100 * ac[1] / pr[1]; sigdigits = C_PR_SIGDIG) ]
    d[2, :] .= [ "true", pr[2], ac[2], round(100 * ac[2] / pr[2]; sigdigits = C_PR_SIGDIG) ]
    d[3, :] .= [ "total", sum(pr), sum(ac), round(100 * sum(ac) / sum(pr); sigdigits = C_PR_SIGDIG) ]
    ## display
    PrettyTables.pretty_table(d, h)
end

## Gelman proposal summary
function tabulate_proposals(results::MultiMCMCResults)
    println("Proposal summary:")
    h = ["Adapted", "Proposed", "Accepted", "Rate"]
    n_iter = length(results.mcmc[1].mc_accepted)
    adp_prd = results.mcmc[1].adapt_period
    ## summarise:
    pr = [length(results.mcmc) * adp_prd, length(results.mcmc) * (n_iter - adp_prd)]
    ac = zeros(Int64, 2)
    for mc in eachindex(results.mcmc)
        ac[1] += sum(results.mcmc[mc].mc_accepted[1:adp_prd])
        ac[2] += sum(results.mcmc[mc].mc_accepted[(adp_prd + 1):n_iter])
    end
    d = Matrix(undef, 3, 4)
    d[1, :] .= [ "false", pr[1], ac[1], round(100 * ac[1] / pr[1]; sigdigits = C_PR_SIGDIG) ]
    d[2, :] .= [ "true", pr[2], ac[2], round(100 * ac[2] / pr[2]; sigdigits = C_PR_SIGDIG) ]
    d[3, :] .= [ "total", sum(pr), sum(ac), round(100 * sum(ac) / sum(pr); sigdigits = C_PR_SIGDIG) ]
    ## display
    PrettyTables.pretty_table(d, h)
end

## results summary
function tabulate_results(results::MCMCResults, proposals = false)
    ## proposals
    proposals && tabulate_proposals(results)
    ## samples
    println("MCMC summary:")
    h = ["θ", "Rμ", "σ", "z"]#, "Iμ"
    d = Matrix(undef, length(results.mean), 4)
    sigma = zeros(length(results.mean))
    for p in eachindex(sigma)
        sigma[p] = sqrt(results.covar[p,p])
    end
    d[:,1] .= 1:length(results.mean)
    d[:,2] .= round.(results.mean; sigdigits = C_PR_SIGDIG)
    d[:,3] .= round.(sigma; sigdigits = C_PR_SIGDIG)
    d[:,4] .= round.(results.geweke[2][1,:]; sigdigits = C_PR_SIGDIG)
    PrettyTables.pretty_table(d, h)
end

## results summary
"""
    tabulate_results

**Parameters**
- `results`     -- the results of a call to `run_multi_chain_analysis`.
- `proposals`   -- display proposal analysis.

Display the results of an MCMC analysis.
"""
function tabulate_results(results::MultiMCMCResults, proposals = false)
    ## proposals
    proposals && tabulate_proposals(results)
    ## samples
    println("Gelman diagnostic:")
    h = ["θ", "Rμ", "σ", "SRE", "SRE95"]#, "Iμ"
    d = Matrix(undef, length(results.mu), 5)
    d[:,1] .= 1:length(results.mu)
    d[:,2] .= round.(results.mu; sigdigits = C_PR_SIGDIG)
    d[:,3] .= round.(results.sdw; sigdigits = C_PR_SIGDIG)
    d[:,4] .= round.(results.sre; sigdigits = C_PR_SIGDIG)
    d[:,5] .= round.(results.sre_ul; sigdigits = C_PR_SIGDIG)
    PrettyTables.pretty_table(d, h)
end
