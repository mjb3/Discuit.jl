# Discuit.jl examples

The following examples provide a flavour of Discuit's core functionality. See the [Discuit.jl manual](@ref) for more detailed instructions.

## SIS model

The following example is based on that published by Pooley et al. in 2015 in the paper that introduces the model based proposal method. We could use `generate_model("SIS", [100,1])` to generate the model but constructing it manually is a helpful exercise for getting to know the package. We start by examining [DiscuitModel](@ref) in the package documentation:

ADD CODE BLOCK TO SEARCH?

Now that we know the necessary parameters for defining a model we can begin by defining a rate function. Note that the correct signature must be used in order for it to be compatible with the package:

```@repl 1
function sis_rf(output::Array{Float64, 1}, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    output[1] = parameters[1] * population[1] * population[2]
    output[2] = parameters[2] * population[2]
end
```

Next we define a simple observation function, again with the correct signature:
```@repl 1
obs_fn(population::Array{Int64, 1}) = population
```

Naturally we choose the same prior distribution as Pooley so that we can compare results. The return type must be Float.

```@repl 1
function weak_prior(parameters::Array{Float64, 1})
    parameters[1] > 0.0 || return 0.0
    parameters[2] > 0.0 || return 0.0
    return 1.0
end
```

Finally, we define an observation likelihood model. Again, we use the same as Pooley, with observation errors normally distributed around the true value with standard deviation `2`:

```@repl 1
function si_gaussian(y::Array{Int64, 1}, population::Array{Int64, 1})
    obs_err = 2
    tmp1 = log(1 / (sqrt(2 * pi) * obs_err))
    tmp2 = 2 * obs_err * obs_err
    obs_diff = y[2] - population[2]
    return tmp1 - ((obs_diff * obs_diff) / tmp2)
end
```
We can now define a model. The three parameters declared inline are the transition matrix; an optional index for the t0 parameter (ignore for now); and the initial condition which represents the state of the population at the origin of each trajectory:

```@repl 1
using Discuit
model = DiscuitModel("SIS", sis_rf, [-1 1; 1 -1], 0, [100, 1], obs_fn, weak_prior, si_gaussian)
```

### Simulation

Although our main goal is to replicate the analysis of Pooley et al. we can also run a simulation using the [gillespie_sim](@ref) function.

```@repl 1
xi = gillespie_sim(model, [0.003,0.1]);
```

We can also visualise the results using the corresponding R package: rDiscuit. ADD LINK

![SIS simulation](https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/sis-sim.png)

### MCMC

Running an MCMC analysis based on a set of observations data is simple. TBC...

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/sis-sim.png" alt="SIS simulation" height="180"/>
```

## Custom MCMC

Some situations...


- link to [Discuit.jl documentation](@ref)
- link to [`set_random_seed(seed::Int64)`](@ref)
