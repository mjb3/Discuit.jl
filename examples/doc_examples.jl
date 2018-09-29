## import the package
using Discuit

### pooley ###
## define model
# rate function
function sis_rf(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    output[1] = parameters[1] * population[1] * population[2]
    output[2] = parameters[2] * population[2]
end
# define obs function (no error)
obs_fn(population::Array{Int64, 1}) = population
# prior
function weak_prior(parameters::Array{Float64, 1})
    parameters[1] > 0.0 || return 0.0
    parameters[2] > 0.0 || return 0.0
    return 1.0
end
# obs model
function si_gaussian(y::Array{Int64, 1}, population::Array{Int64, 1})
    obs_err = 2
    tmp1 = log(1 / (sqrt(2 * pi) * obs_err))
    tmp2 = 2 * obs_err * obs_err
    obs_diff = y[2] - population[2]
    return tmp1 - ((obs_diff * obs_diff) / tmp2)
end
# define model
model = DiscuitModel("SIS", sis_rf, [-1 1; 1 -1], 0, [100, 1], obs_fn, weak_prior, si_gaussian)

## run sim
xi = gillespie_sim(model, [0.003,0.1])

print(length(xi.trajectory))
