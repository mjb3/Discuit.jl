# Discuit.jl

***Fast parameter inference for discrete state space continuous time (DSSCT) models in Julia.***

> Please note that this package is still in development.

## Package features

```@docs
Discuit
```

## Contents

```@contents
```

## Installation

The package can be installed by typing `]` in the REPL to enter the Pkg mode and running:

```
pkg> add https://github.com/mjb3/Discuit.jl
```

## Getting started

The following code initialises a `DiscuitModel` and runs a simulation, storing the results in `x`.

```@repl 1
using Discuit;
set_random_seed(1) # hide
model = generate_model("SIS", [100,1]);
x = gillespie_sim(model, [0.003, 0.1]);
```

We can now run an MCMC analysis using observations data from `x`:

```@repl 1
s = run_met_hastings_mcmc(model, x.observations, [0.0025, 0.12]);
```

## Further usage

More examples can be found in the section [Discuit.jl examples](@ref), including enough code to get up and running with convergence diagnostics and customised models. A more detailed guide to the pre defined models is available in the [Discuit.jl models](@ref) section. Further information regarding the packages other functionality can be found in the [Discuit.jl manual](@ref).
