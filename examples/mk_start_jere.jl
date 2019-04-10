### HMM trade list example

## define a model
function pooley_sis_model()
    ## initial_condition
    ic::Array{Int64, 1} = [100, 1]
    ## rate function
    function sis_rf(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
        output[1] = parameters[1] * population[1] * population[2]
        output[2] = parameters[2] * population[2]
    end
    ## transition matrix
    m_t::Array{Int64, 2} = [-1 1; 1 -1]
    ## define obs function (no error)
    function obs_fn(population::Array{Int64, 1})
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
    function si_gaussian(y::Array{Int64, 1}, population::Array{Int64, 1})
        obs_diff = y[2] - population[2]
        return tmp1 - ((obs_diff * obs_diff) / tmp2)
    end
    ## define model
    comps = "ABCDEFGH"
    return DiscuitModel(comps, ic, sis_rf, m_t, obs_fn, weak_prior, si_gaussian, 0)
end

## chaff
function chaff(a::Array{Int64})
    ndims(a) > 1 && throw("can't handle dim >1") 
    return sum(a)
end

println(chaff([1, 2, 3]))
println(chaff([1 2 3; 1 2 3]))
