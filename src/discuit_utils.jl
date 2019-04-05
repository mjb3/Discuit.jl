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
    print_gelman_results(results::GelmanResults, dpath::String)

**Parameters**
- `results` -- `GelmanResults` variable.
- `dpath`   -- the path of the directory where the results will be saved.

Save the results from a call to `run_gelman_diagnostic` to the directory `dpath`, e.g. "./out/gelman/".
"""
function print_gelman_results(results::GelmanResults, dpath::String)
    # create directory
    # dpath = string("./out/", dname, "/")
    isdir(dpath) || mkpath(dpath)
    # print summary by theta row
    open(string(dpath, "gelman.csv"), "w") do f
        # print headers
        write(f, "parameter,mu,sdw,sre,sre_ll,sre_ul")
        for p in eachindex(results.mu)
            write(f, "\n$p,$(results.mu[p]),$(results.sdw[p]),$(results.sre[p]),$(results.sre_ll[p]),$(results.sre_ul[p])")
        end
    end # end of print summary
    # print chains
    for i in eachindex(results.mcmc)
        print_mcmc_results(results.mcmc[i], string(dpath, "mc", i, "/"))
    end
end

## print MCMC results to file
"""
    print_mcmc_results(mcmc, dpath)

**Parameters**
- `results` -- `MCMCResults` variable.
- `dpath`   -- the path of the directory where the results will be saved.

Save the results from a call to `run_met_hastings_mcmc` or `run_custom_mcmc` to the directory `dpath`, e.g. "./out/mcmc/".
"""
function print_mcmc_results(mcmc::MCMCResults, dpath::String)
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

## proposal summary
function tabulate_proposals(results::GelmanResults)
    println("Proposal summary:")
    h = ["Adapted", "Chain", "Proposed", "Accepted", "Rate"]
    d = Matrix(undef, length(results.mcmc) * 2, 5)
    for mc in eachindex(results.mcmc)
        # burnin
        pr = results.mcmc[mc].adapt_period
        ac = sum(results.mcmc[mc].mc_accepted[1:results.mcmc[mc].adapt_period])
        d[mc, :] .= [ "false", mc, pr, ac, round(100 * ac / pr; sigdigits = C_PR_SIGDIG) ]
    end
    for mc in eachindex(results.mcmc)
        # adapted
        pr = length(results.mcmc[mc].mc_accepted) - results.mcmc[mc].adapt_period
        ac = sum(results.mcmc[mc].mc_accepted[(results.mcmc[mc].adapt_period + 1):length(results.mcmc[mc].mc_accepted)])
        d[(length(results.mcmc) + mc), :] .= [ "true", mc, pr, ac, round(100 * ac / pr; sigdigits = C_PR_SIGDIG) ]
    end
    PrettyTables.pretty_table(d, h)
end

## results summary AS DataFrame?
"""
    tabulate_gelman_results

**Parameters**
- `results`     -- the results of a call to `run_gelman_diagnostic`.
- `proposals`   -- display proposal analysis.

Display the results of a multi chain analysis run using `run_gelman_diagnostic`.
"""
function tabulate_gelman_results(results::GelmanResults, proposals = false)
    ## proposals
    proposals && tabulate_proposals(results)
    ## samples
    println("Gelman diagnostic:")
    h = ["θ", "μ", "σ", "SRE", "SRE95"]
    d = Matrix(undef, length(results.mu), 5)
    d[:,1] .= 1:length(results.mu)
    d[:,2] .= round.(results.mu; sigdigits = C_PR_SIGDIG)
    d[:,3] .= round.(results.sdw; sigdigits = C_PR_SIGDIG)
    d[:,4] .= round.(results.sre; sigdigits = C_PR_SIGDIG)
    d[:,5] .= round.(results.sre_ul; sigdigits = C_PR_SIGDIG)
    PrettyTables.pretty_table(d, h)
end
