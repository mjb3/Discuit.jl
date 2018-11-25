### SIS model from Pooley et al. 2015 ###
# this code defines a SIS model and stores it in the variable 'model'


# rate function
function sis_rf(output::Array{Float64, 1}, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    output[1] = parameters[1] * population[1] * population[2]
    output[2] = parameters[2] * population[2]
end

# transition matrix
t_matrix = [-1 1; 1 -1]

# observation function
obs_fn(population::Array{Int64, 1}) = population

# prior is flat and improper
function weak_prior(parameters::Array{Float64, 1})
    parameters[1] > 0.0 || return 0.0
    parameters[2] > 0.0 || return 0.0
    return 1.0
end

# observation model
function si_gaussian(y::Array{Int64, 1}, population::Array{Int64, 1})
    obs_err = 2
    tmp1 = log(1 / (sqrt(2 * pi) * obs_err))
    tmp2 = 2 * obs_err * obs_err
    obs_diff = y[2] - population[2]
    return tmp1 - ((obs_diff * obs_diff) / tmp2)
end

# initial condition
initial_condition = [100, 1]

# define the model:
model = DiscuitModel("SIS", initial_condition, sis_rf, t_matrix, obs_fn, weak_prior, si_gaussian, 0);
