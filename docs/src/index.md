# Discuit.jl documentation

*Fast parameter inference for discrete state space continuous time models in Julia.*

Discuit: simulation and parameter inference for discrete state space continuous time (DSSCT) models.

!!! note

    Please note that this package is still in development.

## Package Features

    <!-- - Write all your documentation in [Markdown](https://en.wikipedia.org/wiki/Markdown). -->
    - Minimal configuration.
    - Supports Julia `0.6` and `0.7-dev`.
    - Generates tables of contents and docstring indexes.
    - Use `git push` to automatically build and deploy docs from Travis to GitHub Pages.

    The [Discuit.jl documentation](@ref) provides a tutorial explaining how to get started using Documenter.

```@contents
```

## Functions

```@docs
set_random_seed(seed::Int64)
gillespie_sim(model::DiscuitModel, parameters::Array{Float64,1}, tmax::Float64 = 100.0, num_obs::Int64 = 5)
run_met_hastings_mcmc(model::DiscuitModel, obs_data::ObsData, initial_parameters::Array{Float64, 1}, steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
run_custom_mcmc(model::DiscuitModel, obs_data::ObsData, proposal_function::Function, x0::MarkovState, steps::Int64 = 50000, adapt_period::Int64 = 10000, prop_param::Bool = false, ppp::Float64 = 0.3)
```

- link to [Discuit.jl documentation](@ref)
- link to [`set_random_seed(seed::Int64)`](@ref)

## Index

```@index
```
