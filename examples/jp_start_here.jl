### Discuit.jl Pooley (2015) SIS example

### getting started

## install Discuit:
# 1) press `]` to enter Pkg mode
# 2) run: add https://github.com/mjb3/Discuit.jl

## import the package
using Discuit

## define a model
function pooley_sis_model()
    ## initial_condition
    ic::Array{Int32, 1} = [100, 1]
    ## rate function
    function sis_rf(output, parameters::Array{Float64, 1}, population::Array{Int32, 1})
        output[1] = parameters[1] * population[1] * population[2]
        output[2] = parameters[2] * population[2]
    end
    ## transition matrix
    m_t::Array{Int32, 2} = [-1 1; 1 -1]
    ## define obs function (no error)
    function obs_fn(population::Array{Int32, 1})
        return population
    end
    ## prior
    function weak_prior(parameters::Array{Float64, 1})
        parameters[1] > 0.0 || return 0.0
        parameters[2] > 0.0 || return 0.0
        return 1.0
    end
    ## obs model
    obs_err = 2
    # do some preliminary computation
    # - IS THIS the right place to put this?*
    tmp1 = log(1 / (sqrt(2 * pi) * obs_err))
    tmp2 = 2 * obs_err * obs_err
    # define log likelihood function
    function si_gaussian(y::Array{Int32, 1}, population::Array{Int32, 1})
        obs_diff = y[2] - population[2]
        return tmp1 - ((obs_diff * obs_diff) / tmp2)
    end
    ## define model
    # this is a bit of a mess atm, apologies!
    return DiscuitModel("SIS", sis_rf, m_t, 0, ic, obs_fn, weak_prior, si_gaussian)
end

### do some stuff ###
## gillespie example
function sim_example()
    model = pooley_sis_model()
    xi = gillespie_sim(model, [0.003,0.1])
    print_trajectory(model, xi, "./out/sis_sim.csv")
    print_observations(xi.observations, "./out/sis_obs.csv")
end
## run an mcmc analysis and print the results
function mcmc_example()
    # obs data
    obs = Observations([20, 40, 60, 80, 100], [0 18; 0 65; 0 70; 0 66; 0 67])
    model = pooley_sis_model()
    rs = run_met_hastings_mcmc(model, obs, [0.003,0.1])
    print_mcmc_results(rs, "./out/mcmc_example/")
end
## convergence diagnostic
function run_conv_diag()
    obs = Observations([20, 40, 60, 80, 100], [0 18; 0 65; 0 70; 0 66; 0 67])
    model = pooley_sis_model()
    rs = run_gelman_diagnostic(model, obs, [0.0025 0.15; 0.004 0.08; 0.0033 0.1])
    ac = compute_autocorrelation(rs.mcmc)
    print_gelman_results(rs, "./out/gelman_example/")
    print_autocorrelation(ac, "./out/gelman_example/acp.csv")
end

## run the examples
sim_example()
mcmc_example()
run_conv_diag()
