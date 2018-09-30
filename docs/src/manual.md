# Discuit.jl manual

This section contains a directory of `struct`s and `Function`s in the package. See [Discuit.jl examples](@ref) for a broad overview of the package's core functionality.

## Contents

```@contents
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
set_random_seed(seed::Int64)
gillespie_sim(model::DiscuitModel, parameters::Array{Float64,1}, tmax::Float64 = 100.0, num_obs::Int64 = 5)
run_met_hastings_mcmc(model::DiscuitModel, obs_data::Observations, initial_parameters::Array{Float64, 1}, steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
run_custom_mcmc(model::DiscuitModel, obs_data::Observations, proposal_function::Function, x0::MarkovState, steps::Int64 = 50000, adapt_period::Int64 = 10000, prop_param::Bool = false, ppp::Float64 = 0.3)
run_gelman_diagnostic(m_model::DiscuitModel, obs_data::Observations, initial_parameters::Array{Float64, 2}, steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
```

### model helpers

[Discuit.jl](@ref) includes tools for generating components which can help minimise the amount of work required to generate customised [DiscuitModel](@ref)s, including `generate_model(...)` which is used to access a library of pre defined [Discuit.jl models](@ref).

```@docs
generate_generic_obs_function()
generate_weak_prior(n::Int)
generate_gaussian_obs_model(n::Int, σ::AbstractFloat = 2.0)
generate_model(model_name::String, initial_condition::Array{Int64, 1}, obs_error::AbstractFloat = 2.0)
```

### utilities

```@docs
print_trajectory(model::DiscuitModel, sim_results::SimResults, fpath::String)
print_observations(obs_data::Observations, fpath::String)
read_obs_data_from_file(fpath::String)
print_mcmc_results(mcmc::MCMCResults, dpath::String)
print_gelman_results(results::GelmanResults, dpath::String)
```

### custom MCMC

TBC...

## Index

```@index
```

```@raw html
<html>
<head>
<script type="text/javascript" src="https://cdn.rawgit.com/pcooksey/bibtex-js/ef59e62c/src/bibtex_js.js"></script>
</head>
<body>
<textarea id="bibtex_input" style="display:none;">
@book{book1,
  author = "Donald Knuth",
  title = "Concrete Mathematics"
}
</textarea>
<div id="bibtex_display"></div>
</body>
</html>
```
