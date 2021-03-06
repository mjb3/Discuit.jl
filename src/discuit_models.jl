### pre defined models for Discuit.jl ###
## function generating functions
# weak prior
"""
    generate_weak_prior(n)

**Parameters**

- `n`   -- the number of parameters in the model.
- `b`   -- the upper bound of the [Uniform] distribution.

# Examples

    generate_weak_prior(1)

Generate a "weak" prior distribution, Uniform multivariate ~ U(0, max) for dim = n,  where `n` is the number of parameters in the model.
"""
function generate_weak_prior(n::Int, b::Float64 = 1.0)
    # function priorr_density(parameters::Array{Float64, 1})
    #     for i in eachindex(parameters)
    #         parameters[i] < 0.0 && (return 0.0)
    #     end
    #     return 1.0
    # end
    # return priorr_density
    return Distributions.Product(Distributions.Uniform.(zeros(n), b))
end
# gaussian observation likelihood model
"""
    generate_gaussian_obs_model(n, σ = 2.0)

**Parameters**
- `n`   -- the number of discrete states in the model.
- `σ`   -- observation error.

test latex eqn:

```math
\frac{n!}{k!(n - k)!} = \binom{n}{k}
```

# Examples

    p = generate_gaussian_obs_model(2)

Generate a Gaussian observation model for a model with `n` states. Optionally specify observation error `σ`.
"""
function generate_gaussian_obs_model(n::Int, σ::AbstractFloat = 2.0)
    # do some preliminary computation
    tmp1 = log(1 / (sqrt(2 * pi) * σ))
    tmp2 = 2 * σ * σ
    if n == 2
        function gom1(y::Array{Int64, 1}, population::Array{Int64, 1})
            return tmp1 - ( ( (y[2] - population[2]) ^ 2) / tmp2 )
        end
        return gom1
    else
        # SHOULD I VECTORISE THIS? *********
        function gom2(y::Array{Int64, 1}, population::Array{Int64, 1})
            output = 0.0
            for i in 2:n
                output += tmp1 - ( ( (y[i] - population[i]) ^ 2) / tmp2 )
            end
            return output
        end
        return gom2
    end
end

## rate functions
# SI
function si_rf(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    output[1] = parameters[1] * population[1] * population[2]
end
# SIR; SIS
function sir_rf(output::Array{Float64, 1}, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    output[1] = parameters[1] * population[1] * population[2]
    output[2] = parameters[2] * population[2]
end
# SEI
function sei_rf(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    output[1] = parameters[1] * population[1] * population[3]
    output[2] = parameters[2] * population[2]
end
# SEIR
function seir_rf(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    output[1] = parameters[1] * population[1] * population[3]
    output[2] = parameters[2] * population[2]
    output[3] = parameters[3] * population[3]
end
# Lotka-Volterra
function lotka_rf(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    # prey; predator reproduction; predator death
    output[1] = parameters[1] * population[2]
    output[2] = parameters[2] * population[1] * population[2]
    output[3] = parameters[3] * population[1]
end

## generic observation function (no error)
"""
    generate_generic_obs_function()

Generates a simple observation function for use in a `DiscuitModel`. Not very realistic...
"""
function generate_generic_obs_function()
    obs_fn(population::Array{Int64, 1}) = population
    return obs_fn
end
"""
    generate_model(model_name, initial_condition, σ = 2.0)

**Parameters**
- `model_name`          -- the model, e.g. "SI"; "SIR"; "SEIR"; etc
- `initial_condition`   -- initial condition.
- `σ`                   -- observation error.

`model_name` **options**
- `"SI"`
- `"SIR"`
- `"SIS"`
- `"SEI"`
- `"SEIR"`
- `"SEIS"`
- `"SEIRS"`
- `"PREDPREY"`
- `"ROSSMAC"`

# Examples

    generate_model("SIS", [100,1])

Generates a `DiscuitModel`. Optionally specify observation error `σ`.
"""
function generate_model(model_name::String, initial_condition::Array{Int64, 1}, obs_error::AbstractFloat = 2.0)
    # force to upper case here ***
    obs_model = generate_gaussian_obs_model(length(initial_condition), obs_error)
    if model_name == "SI"
        rate_fn = si_rf
        m_transition = [-1 1;]
    elseif model_name == "SIR"
        rate_fn = sir_rf
        m_transition = [-1 1 0; 0 -1 1]
    elseif model_name == "SIS"
        rate_fn = sir_rf
        m_transition = [-1 1; 1 -1]
    elseif model_name == "SEI"
        rate_fn = sei_rf
        m_transition = [-1 1 0; 0 -1 1]
    elseif model_name == "SEIR"
        rate_fn = seir_rf
        m_transition = [-1 1 0 0; 0 -1 1 0; 0 0 -1 1]
    elseif model_name == "LOTKA"
        model_name = "PN"
        rate_fn = lotka_rf
        m_transition = [0 1; 1 -1; -1 0]
    else
        println("ERROR: ", model_name, " not recognised.")
        # handle this better? ***
        return 0
    end
    # NEED TO ADD MORE MODELS ******************
    prior = generate_weak_prior(size(m_transition, 1))
    return DiscuitModel(model_name, initial_condition, rate_fn, m_transition, generate_generic_obs_function(), prior, obs_model, 0)
end
