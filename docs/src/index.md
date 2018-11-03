# Discuit.jl

***Fast parameter inference for discrete state space continuous time (DSSCT) models in Julia.***

> Please note that this package is still in development.

Discuit is a package for Bayesian inference in Discrete state space continuous time (DSSCT) models. DSSCT models, sometimes referred to as compartmental models, are used to represent systems where individuals are assumed, usually as a simplifying abstraction, to move between discrete states. They are a well established research tool in fields including physics, chemistry and ecology. Their use in scientific research typically involves the challenge of comparing ‘noisy’ experimental or observational data to the unobserved (i.e. latent) underlying processes described by the model. Bayes’ theorem provides a convenient framework for dealing with uncertainty in the data. The Bayesian framework also provide a natural, methodologically consistent way of incorporating existing scientific knowledge of the system; the prior distribution of the model parameters

The augmented data Markov chain Monte Carlo (MCMC) methods implemented in Discuit work by introducing a latent variable $\xi$ which represents the sequence of events in a single realisation of the model:

$\pi(\theta|y) = \pi(y|\xi) \pi(\xi|\theta) \pi(\theta)$

See [Introduction to MCMC](@ref) for a basic introduction to MCMC and [Introduction to Monte Carlo methods](@ref) for an overview of random sampling generally. Two algorithms for making proposals to the augmented data space are shipped with the package, with user defined implementations made possible via an alternative [custom MCMC](@ref) framework. Automated tools for analysis and convergence diagnostics include autocorrelation, the Geweke test of stationarity and the Gelman-Rubin diagnostic for multiple Markov chains (a convenient way to run analyses where more than one processor thread is available for use). [Simulation](@ref) via the Gillespie direct method provides a source of simulated observations data for evaluation and validation of the core inference functionality.

See the [Discuit.jl models](@ref) section for an introduction to the aforementioned functionality and the [Discuit.jl manual](@ref) for a description of data types and functions. See the [Discuit in R](https://mjb3.github.io/Discuit/) package documentation for a description of the equivalent functionality in that package.

## Package features

```@docs
Discuit
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
set_random_seed(1); # hide
model = generate_model("SIS", [100,1]);
x = gillespie_sim(model, [0.003, 0.1]);
```

We can now run an MCMC analysis using observations data from `x`:

```@repl 1
s = run_met_hastings_mcmc(model, x.observations, [0.0025, 0.12]);
```

## Further usage

More examples can be found in the section [Discuit.jl examples](@ref), including enough code to get up and running with convergence diagnostics and customised models. A more detailed guide to the pre defined models is available in the [Discuit.jl models](@ref) section. Further information regarding the packages other functionality can be found in the [Discuit.jl manual](@ref).

## Tutorials

* [Introduction to Monte Carlo methods](https://mjb3.github.io/Discuit/articles/monte_carlo_intro/monte_carlo_intro.html): a beginner's guide in Python.
* A basic [Introduction to MCMC](https://mjb3.github.io/Discuit/articles/mcmc_intro/mcmc_intro.html) methods in Python.
* [Discuit.jl examples](@ref): an introduction to MCMC and simulation in [Discuit.jl](@ref) for Julia.
* See the [Discuit for R package documentation](https://mjb3.github.io/Discuit/articles/examples.html) for R examples.
