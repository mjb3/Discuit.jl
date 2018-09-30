### Discuit custom data structures ###

## event instance
struct Event
    time::Float64
    event_type::Int16
end
## observations data
"""
    Observations

Stores one column vector of observation times and one or more column vectors of observation integer values.

**Fields**
- `time`    -- observation times.
- `val`     -- observation values.
"""
struct Observations
    time::Array{Float64, 1}
    val::Array{Int64, 2}
end
## Discuit model
"""
    DiscuitModel

**Fields**
- `model_name`          -- string, e,g, `"SIR"`.
- `rate_function`       -- event rate function.
- `m_transition`        -- transition matrix.
- `t0_index`            -- index of the parameter that represents the initial time. `0` if fixed at `0.0`.
- `initial_condition`   -- initial condition
- `obs_function`        -- observation function.
- `prior`               -- prior density function.
- `obs_model`           -- observation likelihood model.

A `mutable struct` which represents a DSSCT model for use with [Discuit.jl](@ref) functions.
"""
mutable struct DiscuitModel{RFT<:Function, OFT<:Function, PDT<:Function, OMT<:Function}
    # model name
    model_name::String
    # event rate function
    rate_function::RFT
    # transition matrix
    m_transition::Array{Int64, 2}
    # t0 index (0 ~ fixed at 0.0)
    t0_index::Int64
    # initial condition
    initial_condition::Array{Int64, 1}
    # observation function (for sim TBA)
    obs_function::OFT
    # prior density function
    prior::PDT
    # observation model (log likelihood)
    obs_model::OMT
end
## JIT private models
function get_private_model(model::DiscuitModel, obs_data::Observations)
    return PrivateDiscuitModel(model.rate_function, model.m_transition, model.initial_condition, model.t0_index, model.obs_function, model.prior, model.obs_model, obs_data)
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
    obs_function::OFT
    # prior density function
    prior::PDT
    # observation model (log likelihood)
    obs_model::OMT
    # obs data
    obs_data::Observations
end
# results of gillespie sim
"""
    SimResults

**Fields**
- `trajectory`      -- array of type `Event`.
- `observations`    -- variable of type `Observations`.

The results of a simulation.
"""
struct SimResults
    trajectory::Array{Event}
    observations::Observations
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
    trajectory::Array{Event, 1}
    log_like::Float64
    prop_like::Float64
    # TO BE REMOVED?
    prop_type::Int64
end
# results of an McMC analysis
"""
    MCMCResults

**Fields**
- `samples`     -- two dimensional array of samples.
- `mc_accepted` -- proposal accepted
- `mean`        -- sample mean.
- `covar`       -- parameter covariance matrix.

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
end
# NEED TO ADD A VARIANCE MEASURE ***
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
