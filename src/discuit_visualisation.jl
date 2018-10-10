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
function plot_parameter_trace(mcmc::MCMCResults, parameter::Int64)
    x = mcmc.adapt_period:size(mcmc.samples, 1)
    PyPlot.plot(x, mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), parameter])
    PyPlot.xlabel("sample")
    PyPlot.ylabel(string("\$\\theta_", parameter, "\$"))
end
## marginal
function plot_parameter_marginal(mcmc::MCMCResults, parameter::Int64)
    x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), parameter]
    PyPlot.plt[:hist](x, 30)
    # PyPlot.title("TITLE")
    PyPlot.xlabel(string("\$\\theta_", parameter, "\$"))
    PyPlot.ylabel("density")
end
## heatmap
function plot_parameter_heatmap(mcmc::MCMCResults, x_parameter::Int64, y_parameter::Int64)
    x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), x_parameter]
    y = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), y_parameter]
    PyPlot.hexbin(x, y)
    PyPlot.xlabel(string("\$\\theta_", x_parameter, "\$"))
    PyPlot.ylabel(string("\$\\theta_", y_parameter, "\$"))
end
