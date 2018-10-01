# Discuit.jl manual

This section contains a directory of `struct`s and `Function`s in the package. See [Discuit.jl examples](@ref) for a broad overview of the package's core functionality.

## Contents

```@contents
Pages = ["manual.md"]
Depth = 2
```

## Custom structs

```@docs
DiscuitModel
SimResults
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
run_custom_mcmc
run_gelman_diagnostic
compute_autocorrelation
```

### model helpers

[Discuit.jl](@ref) includes tools for generating components which can help minimise the amount of work required to generate customised [DiscuitModel](@ref)s, including `generate_model(...)` which is used to access a library of pre defined [Discuit.jl models](@ref).

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
read_obs_data_from_file
print_mcmc_results
print_gelman_results
```

### custom MCMC

TBC...

## Index

```@index
```

## References

TBA
