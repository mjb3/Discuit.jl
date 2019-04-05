## https://github.com/Evizero/UnicodePlots.jl
import UnicodePlots

## trajectory
"""
    plot_trajectory(x)

**Parameters**
- `x`       -- `SimResults`, i.e. from a call to `gillespie_sim`.

Plot the trajectory of a a DGA simulation on `model` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
"""
function plot_trajectory(x::SimResults)
    p = UnicodePlots.lineplot(x.trajectory.time, x.population[:,1], title = string(x.model_name, " simulation"), name = string(x.model_name[1]))
    for i in 2:size(x.population, 2)
        UnicodePlots.lineplot!(p, x.trajectory.time, x.population[:,i], name = string(x.model_name[i]))
    end
    UnicodePlots.xlabel!(p, "time")
    UnicodePlots.ylabel!(p, "population")
    return p
end

## traceplot
# single
"""
    plot_parameter_trace(mcmc, parameter)

**Parameters**
- `mcmc`        -- `MCMCResults`, e.g. from a call to `run_met_hastings_mcmc`.
- `parameter`   -- the index of the model parameter to be plotted.

Trace plot of samples from an MCMC analysis for a given model `parameter` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
"""
function plot_parameter_trace(mcmc::MCMCResults, parameter::Int64)
    x = 1:size(mcmc.samples, 1)
    p = UnicodePlots.lineplot(x, mcmc.samples[:, parameter], title = string("θ", Char(8320 + parameter), " traceplot."))
    UnicodePlots.xlabel!(p, "sample")
    UnicodePlots.ylabel!(p, string("θ", Char(8320 + parameter)))
    return p
end
# multiple
"""
    plot_parameter_trace(mcmc, parameter)

**Parameters**
- `mcmc`        -- array of `MCMCResults`, e.g. from a call to `run_gelman_diagnostic`.
- `parameter`   -- the index of the model parameter to be plotted.

Trace plot of samples from `n` MCMC analyses for a given model `parameter` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
"""
function plot_parameter_trace(mcmc::Array{MCMCResults, 1}, parameter::Int64)
    x = 1:size(mcmc[1].samples, 1)
    p = UnicodePlots.lineplot(x, mcmc[1].samples[:, parameter], title = string("θ", Char(8320 + parameter), " traceplot."))
    for i in 2:length(mcmc)
        UnicodePlots.lineplot!(p, mcmc[i].samples[:, parameter])
    end
    UnicodePlots.xlabel!(p, "sample")
    UnicodePlots.ylabel!(p, string("θ", Char(8320 + parameter)))
    return p
end

## marginal
"""
    plot_parameter_marginal(mcmc, parameter)

**Parameters**
- `mcmc`        -- `MCMCResults`, e.g. from a call to `run_met_hastings_mcmc`.
- `parameter`   -- the index of the model parameter to be plotted.

Plot the marginal distribution of samples from an MCMC analysis for a given model `parameter` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
"""
function plot_parameter_marginal(mcmc::MCMCResults, parameter::Int64)
    x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), parameter]
    p = UnicodePlots.histogram(x, bins = 20)
    UnicodePlots.ylabel!(p, string("θ", Char(8320 + parameter)))
    UnicodePlots.xlabel!(p, "samples")
    return p
end

## heatmap
"""
    plot_parameter_heatmap(mcmc, x_parameter, y_parameter)

**Parameters**
- `mcmc`        -- `MCMCResults`, e.g. from a call to `run_met_hastings_mcmc`.
- `x_parameter`   -- the index of the model parameter to be plotted on the x axis.
- `y_parameter`   -- the index of the model parameter to be plotted on the y axis.

Plot the marginal distribution of samples from an MCMC analysis for two model parameters using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
"""
function plot_parameter_heatmap(mcmc::MCMCResults, x_parameter::Int64, y_parameter::Int64)
    x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), x_parameter]
    y = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), y_parameter]
    p = UnicodePlots.densityplot(x, y, color = :red)
    UnicodePlots.xlabel!(p, string("θ", Char(8320 + x_parameter)))
    UnicodePlots.ylabel!(p, string("θ", Char(8320 + y_parameter)))
    return p
end

## geweke
"""
    plot_geweke_series(mcmc)

**Parameters**
- `mcmc`        -- `MCMCResults`, e.g. from a call to `run_met_hastings_mcmc`.

Plot the Geweke series...
"""
function plot_geweke_series(mcmc::MCMCResults)
    x = mcmc.geweke[1]
    p = UnicodePlots.scatterplot(x, mcmc.geweke[2][:,1])
    for i in 2:size(mcmc.geweke[2], 2)
        UnicodePlots.scatterplot!(p, x, mcmc.geweke[2][:,i])
    end
    UnicodePlots.lineplot!(p, -2.0, 0.0, color = :yellow)
    UnicodePlots.lineplot!(p, 2.0, 0.0, color = :yellow)
    UnicodePlots.ylabel!(p, "z")
    return p
end

## autocorrelation R
"""
    plot_autocorrelation(autocorrelation)

**Parameters**
- `autocorrelation`     -- The results from a call to `compute_autocorrelation`.

Plot autocorrelation for an MCMC analysis.
"""
function plot_autocorrelation(autocorrelation::AutocorrelationResults)
    # build y
    for i in eachindex(autocorrelation.autocorrelation)
        autocorrelation.autocorrelation[i] = max(autocorrelation.autocorrelation[i], 0)
    end
    # plot
    p = UnicodePlots.lineplot(autocorrelation.lag, autocorrelation.autocorrelation[:,1], title = string("θ autocorrelation"))
    for i in 2:size(autocorrelation.autocorrelation, 2)
        UnicodePlots.lineplot!(p, autocorrelation.lag, autocorrelation.autocorrelation[:,i])
    end
    UnicodePlots.xlabel!(p, "lag")
    UnicodePlots.ylabel!(p, "R")
    return p
end
# function plot_autocorrelation(autocorrelation::Array{Float64, 2})
#     # build x
#     x = zeros(size(autocorrelation, 1))
#     for i in eachindex(x)
#         x[i] = (i - 1) * AC_LAG_INT
#     end
#     # build y
#     for i in eachindex(autocorrelation)
#         autocorrelation[i] = max(autocorrelation[i], 0)
#     end
#     # plot
#     p = UnicodePlots.lineplot(x, autocorrelation[:,1], title = string("θ autocorrelation"))
#     for i in 2:size(autocorrelation, 2)
#         UnicodePlots.lineplot!(p, x, autocorrelation[:,i])
#     end
#     UnicodePlots.xlabel!(p, "lag")
#     UnicodePlots.ylabel!(p, "R")
#     return p
# end

# single
# """
#     plot_autocorrelation(mcmc)
#
# **Parameters**
# - `mcmc`        -- `MCMCResults`, e.g. from a call to `run_met_hastings_mcmc`.
#
# Plot autocorrelation R using ADD PYPLOT LINK.
# """
# function plot_autocorrelation(mcmc::MCMCResults)
#
# end
