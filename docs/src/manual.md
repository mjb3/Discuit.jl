# Discuit.jl manual

See [Discuit.jl examples](@ref) for a brief introduction to the package's core functionality.

```@contents
Pages = ["manual.md"]
Depth = 3
```

## Types

```@docs
DiscuitModel
SimResults
Trajectory
Observations
MCMCResults
GelmanResults
```

## Functions

This section is organised in three parts:
- the main package [core functionality](@ref) for working with standard Discuit models
- [utilities](@ref), for loading to and from file
- [custom MCMC](@ref), for running custom algorithms

### core functionality

```@docs
set_random_seed
gillespie_sim
run_met_hastings_mcmc
run_gelman_diagnostic
compute_autocorrelation
```

### model helpers

[Discuit.jl](@ref) includes tools for generating components which can help minimise the amount of work required to generate customised `DiscuitModel`s, including `generate_model(...)` which is used to access a library of pre defined [Discuit.jl models](@ref).

```@docs
generate_generic_obs_function
generate_weak_prior(n::Int)
generate_gaussian_obs_model(n::Int, σ::AbstractFloat = 2.0)
generate_model
```

### utilities

```@docs
print_trajectory
print_observations
get_observations_from_file
get_observations_from_array
print_mcmc_results
print_gelman_results
print_autocorrelation
```

### visualisation

```@docs
plot_trajectory
plot_parameter_trace
plot_parameter_marginal
plot_parameter_heatmap
plot_geweke_series
plot_autocorrelation
```

### custom MCMC

NEED to add: generate_custom_x0

```@docs
run_custom_mcmc
run_custom_mcmc_gelman_diagnostic
```

## Index

```@index
```

## References

TBA
