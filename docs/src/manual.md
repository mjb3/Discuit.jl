# Discuit.jl manual

## Custom structs

```@docs
DiscuitModel
MCMCResults
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
run_met_hastings_mcmc(model::DiscuitModel, obs_data::ObsData, initial_parameters::Array{Float64, 1}, steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
run_custom_mcmc(model::DiscuitModel, obs_data::ObsData, proposal_function::Function, x0::MarkovState, steps::Int64 = 50000, adapt_period::Int64 = 10000, prop_param::Bool = false, ppp::Float64 = 0.3)
run_gelman_diagnostic(m_model::DiscuitModel, obs_data::ObsData, initial_parameters::Array{Float64, 2}, steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
```

### utilities
```@docs
print_trajectory(model::DiscuitModel, sim_results::SimResults, fpath::String)
print_observations(obs_data::ObsData, fpath::String)
read_obs_data_from_file(fpath::String)
print_mcmc_results(mcmc::McMCResults, dpath::String)
print_gelman_results(results::GelmanResults, dpath::String)
```

### custom MCMC
TBC...
