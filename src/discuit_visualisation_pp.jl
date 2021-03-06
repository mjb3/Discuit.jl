## PYPLOT
import PyPlot

## trajectory
"""
    plot_trajectory(x)

**Parameters**
- `x`       -- `SimResults`, i.e. from a call to `gillespie_sim`.

Plot the trajectory of a a DGA simulation on `model` using ADD PYPLOT LINK.
"""
function plot_trajectory(x::SimResults)
    PyPlot.plot(x.trajectory.time, x.population)
    PyPlot.xlabel("time")
    PyPlot.ylabel("population")
end

## traceplot
# single
"""
    plot_parameter_trace(mcmc, parameter)

**Parameters**
- `mcmc`        -- `MCMCResults`, e.g. from a call to `run_single_chain_analysis`.
- `parameter`   -- the index of the model parameter to be plotted.

Trace plot of samples from an MCMC analysis for a given model `parameter` using ADD PYPLOT LINK.
"""
function plot_parameter_trace(mcmc::MCMCResults, parameter::Int64)
    x = mcmc.adapt_period:size(mcmc.samples, 1)
    PyPlot.plot(x, mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), parameter])
    PyPlot.xlabel("sample")
    PyPlot.ylabel(string("\$\\theta_", parameter, "\$"))
end
# multiple
"""
    plot_parameter_trace(mcmc, parameter)

**Parameters**
- `mcmc`        -- array of `MCMCResults`, e.g. from a call to `run_multi_chain_analysis`.
- `parameter`   -- the index of the model parameter to be plotted.

Trace plot of samples from `n` MCMC analyses for a given model `parameter` using ADD PYPLOT LINK.
"""
function plot_parameter_trace(mcmc::Array{MCMCResults, 1}, parameter::Int64)
    x = mcmc[1].adapt_period:size(mcmc[1].samples, 1)
    for i in eachindex(mcmc)
        PyPlot.plot(x, mcmc[i].samples[mcmc[i].adapt_period:size(mcmc[i].samples, 1), parameter])
    end
    PyPlot.xlabel("sample")
    PyPlot.ylabel(string("\$\\theta_", parameter, "\$"))
end
## marginal
"""
    plot_parameter_marginal(mcmc, parameter)

**Parameters**
- `mcmc`        -- `MCMCResults`, e.g. from a call to `run_single_chain_analysis`.
- `parameter`   -- the index of the model parameter to be plotted.

Plot the marginal distribution of samples from an MCMC analysis for a given model `parameter` using ADD PYPLOT LINK.
"""
function plot_parameter_marginal(mcmc::MCMCResults, parameter::Int64)
    x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), parameter]
    PyPlot.plt[:hist](x, 30)
    PyPlot.xlabel(string("\$\\theta_", parameter, "\$"))
    PyPlot.ylabel("density")
end
## heatmap
"""
    plot_parameter_heatmap(mcmc, x_parameter, y_parameter)

**Parameters**
- `mcmc`        -- `MCMCResults`, e.g. from a call to `run_single_chain_analysis`.
- `x_parameter`   -- the index of the model parameter to be plotted on the x axis.
- `y_parameter`   -- the index of the model parameter to be plotted on the y axis.

Plot the marginal distribution of samples from an MCMC analysis for two model parameters using ADD PYPLOT LINK.
"""
function plot_parameter_heatmap(mcmc::MCMCResults, x_parameter::Int64, y_parameter::Int64)
    x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), x_parameter]
    y = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), y_parameter]
    PyPlot.hexbin(x, y)
    PyPlot.xlabel(string("\$\\theta_", x_parameter, "\$"))
    PyPlot.ylabel(string("\$\\theta_", y_parameter, "\$"))
end
## autocorrelation R
# single
# """
#     plot_autocorrelation(mcmc)
#
# **Parameters**
# - `mcmc`        -- `MCMCResults`, e.g. from a call to `run_single_chain_analysis`.
#
# Plot autocorrelation R using ADD PYPLOT LINK.
# """
# function plot_autocorrelation(mcmc::MCMCResults)
#
# end
