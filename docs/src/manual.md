# Discuit.jl manual

## Custom structs

## Functions

```@docs
set_random_seed(seed::Int64)
gillespie_sim(model::DiscuitModel, parameters::Array{Float64,1}, tmax::Float64 = 100.0, num_obs::Int64 = 5)
run_met_hastings_mcmc(model::DiscuitModel, obs_data::ObsData, initial_parameters::Array{Float64, 1}, steps::Int64 = 50000, adapt_period::Int64 = 10000, mbp::Bool = true, ppp::Float64 = 0.3)
run_custom_mcmc(model::DiscuitModel, obs_data::ObsData, proposal_function::Function, x0::MarkovState, steps::Int64 = 50000, adapt_period::Int64 = 10000, prop_param::Bool = false, ppp::Float64 = 0.3)
```
