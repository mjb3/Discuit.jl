### Discuit custom data structures ###

## event instance
# struct Event
#     time::Float64
#     event_type::Int16
# end
# trajectory
"""
    Trajectory

**Fields**
- `time`        -- event times.
- `event_type`  -- event type, index of `model.rate_function`.

A single realisation of the model.
"""
struct Trajectory
    time::Array{Float64, 1}
    event_type::Array{Int64, 1}
end
# observations data
"""
    Observations

**Fields**
- `time`    -- observation times.
- `val`     -- observation values.

# Examples

    # pooley dataset
    pooley = Observations([20, 40, 60, 80, 100], [0 18; 0 65; 0 70; 0 66; 0 67])

Stores one column vector of observation times and one or more column vectors of observation integer values.
"""
struct Observations
    time::Array{Float64, 1}
    val::Array{Int64, 2}
end
## results of gillespie sim
"""
    SimResults

**Fields**
- `trajectory`      -- a variable of type `Trajectory`.
- `population`      -- gives the state of the system after each event.
- `observations`    -- variable of type `Observations`.

The results of a simulation.
"""
struct SimResults
    model_name::String
    trajectory::Trajectory
    population::Array{Int64, 2}
    observations::Observations
end

## Discuit model
"""
    DiscuitModel

**Fields**
- `model_name`          -- string, e,g, `"SIR"`.
- `initial_condition`   -- initial condition
- `rate_function`       -- event rate function.
- `m_transition`        -- transition matrix.
- `observation_function -- observation function, use this to add 'noise' to simulated observations.
- `prior_density`       -- prior density function.
- `observation_model`   -- observation model likelihood function.
- `t0_index`            -- index of the parameter that represents the initial time. `0` if fixed at `0.0`.

A `mutable struct` which represents a DSSCT model (see [Discuit.jl models](@ref) for further details).
"""
mutable struct DiscuitModel
    # model name
    model_name::String
    # initial condition
    initial_condition::Array{Int64, 1}
    # event rate function
    rate_function::Function
    # transition matrix
    m_transition::Array{Int64, 2}
    # observation function (for sim TBA)
    observation_function::Function
    # prior density function
    prior_density::Function
    # observation model (log likelihood)
    observation_model::Function
    # t0 index (0 ~ fixed at 0.0)
    t0_index::Int64
end
## JIT private models
function get_private_model(model::DiscuitModel, obs_data::Observations)
    return PrivateDiscuitModel(model.rate_function, model.m_transition, model.initial_condition, model.t0_index, model.observation_function, model.prior_density, model.observation_model, obs_data)
end
# latent observation model
struct PrivateDiscuitModel{RFT<:Function, OFT<:Function, PDT<:Function, OMT<:Function}
    # event rate function
    rate_function::RFT
    # transition matrix
    m_transition::Array{Int64, 2}
    # initial condition
    initial_condition::Array{Int64, 1}
    # t0 index (0 ~ fixed at 0.0)
    t0_index::Int64
    # observation function (for sim TBA)
    observation_function::OFT
    # prior density function
    prior_density::PDT
    # observation model (log likelihood)
    observation_model::OMT
    # obs data
    obs_data::Observations
end
## generic proposal data structures
# parameter
struct ParameterProposal
    value::Array{Float64, 1}
    prior::Float64
end
# state
struct MarkovState
    parameters::ParameterProposal
    trajectory::Trajectory
    log_like::Float64
    prop_like::Float64
    # TO BE REMOVED?
    prop_type::Int64
end

## results of an McMC analysis
"""
    MCMCResults

**Fields**
- `samples`         -- two dimensional array of samples.
- `mc_accepted`     -- proposal accepted
- `mean`            -- sample mean.
- `covar`           -- parameter covariance matrix.
- `proposal_alg`    -- proposal algorithm.
- `num_obs`         -- number of observations.
- `adapt_period`    -- adaptation (i.e. 'burn in') period.
- `geweke`          -- Geweke convergence diagnostic results (`Tuple` of `Array`s).

The results of an MCMC analysis including samples; mean; covariance matrix; adaptation period; and results of the Geweke test of stationarity.
"""
struct MCMCResults
    samples::Array{Float64, 2}
    mc_accepted::Array{Float64, 1}
    mean::Array{Float64, 1}
    covar::Array{Float64, 2}
    proposal_alg::String
    num_obs::Int64
    adapt_period::Int64
    geweke::Tuple{Array{Int64, 1},Array{Float64, 2}}
    # TO BE REMOVED (DEBUG)
    mcf::Array{Float64, 2}  # TO BE REMOVED
    mc_log_like::Array{Float64, 1}  # TO BE REMOVED
    xn::MarkovState
    prop_type::Array{Int64, 1}
    ll_g::Array{Float64, 1}
    mh_prob::Array{Float64, 1}
    mc_time::Array{Float64, 1}  # TO BE REMOVED
end

## autocorrelation results
"""
    AutocorrelationResults

**Fields**
- `lag`             -- autocorrelation lag.
- `autocorrelation` -- autocorrelation statistics.

Results of a call to `compute_autocorrelation`.
"""
struct AutocorrelationResults
    lag::Array{Int64,1}
    autocorrelation::Array{Float64,2}
end

## Gelman Rubin test results
"""
    GelmanResults

**Fields**
- `mu`      -- between chain sample mean.
- `sre`     -- scale reduction factor estimate.
- `sre_ll`  -- scale reduction factor lower confidence interval.
- `sre_ul`  -- scale reduction factor upper confidence interval.
- `mcmc`    -- array of `MCMCResults`

Results of a Gelman Rubin convergence diagnostic including n `MCMCResults` variables; `mu`; and the scale reduction factor estimates (`sre`).
"""
struct GelmanResults
    mu::Array{Float64, 1}
    sre::Array{Float64, 1}
    sre_ll::Array{Float64, 1}
    sre_ul::Array{Float64, 1}
    mcmc::Array{MCMCResults,1}
end
