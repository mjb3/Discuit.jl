# Discuit.jl examples

The following examples provide a flavour of Discuit's core functionality. See the [Discuit.jl manual](@ref) for more detailed instructions.

## SIS model

The following example is based on that published by Pooley et al. in 2015 in the paper that introduces the model based proposal method. To recreate Pooley's analysis we must first define the a [DiscuitModel](@ref):

```@repl
a = 1
b = 2
a + b
```

```@repl
a = 1
b = 2
a + b
```

### Simulation

![SIS simulation](https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/sis-sim.png)

### MCMC

Running an MCMC analysis based on a set of observations data is simple. TBC...

## Custom MCMC

Some situations...


- link to [Discuit.jl documentation](@ref)
- link to [`set_random_seed(seed::Int64)`](@ref)
