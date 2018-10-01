### pre defined models for Discuit.jl ###
## function generating functions
# weak prior
"""
    generate_weak_prior(n)

**Parameters**

- `n`   -- the number of parameters in the model.

# Examples

    generate_weak_prior(1)

Generate a "weak" prior density function, where `n` is the number of parameters in the model.
"""
function generate_weak_prior(n::Int)
    function prior_density(parameters::Array{Float64, 1})
        for i in eachindex(parameters)
            parameters[i] < 0.0 && (return 0.0)
        end
        return 1.0
    end
    return prior_density
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

    p = generate_weak_prior(1)

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
    output[2] = parameters[2] * population[2]
end
# SIR; SIS
function sir_rf(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
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
## generic observation function (no error)
"""
    generate_generic_obs_function()

# Examples

    generate_generic_obs_function()

Generates a simple observation function for use in a [DiscuitModel](@ref). Not very realistic...
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

Generates a [DiscuitModel](@ref). Optionally specify observation error `σ`.
"""
function generate_model(model_name::String, initial_condition::Array{Int64, 1}, obs_error::AbstractFloat = 2.0)
    # force to upper case here? ***
    obs_model = generate_gaussian_obs_model(length(initial_condition), obs_error)
    if model_name == "SI"
        m_transition = [-1 1;]
    elseif model_name == "SIR"
        m_transition = [-1 1 0; 0 -1 1]
    elseif model_name == "SIS"
        m_transition = [-1 1; 1 -1]
    elseif model_name == "SEI"
        m_transition = [-1 1 0; 0 -1 1]
    elseif model_name == "SEIR"
        m_transition = [-1 1 0 0; 0 -1 1 0; 0 0 -1 1]
    else
        print("\nERROR: ", model_name, " not recognised.")
        # handle this better? ***
        return 0
    end
    # NEED TO ADD MORE MODELS ******************
    prior = generate_weak_prior(size(m_transition, 1))
    return DiscuitModel(model_name, sir_rf, m_transition, 0, initial_condition, generate_generic_obs_function(), prior, obs_model)
end
