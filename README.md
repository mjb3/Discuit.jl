# Discuit.jl
Fast parameter inference for discrete state space continuous time (DSSCT) models in Julia

| **Build Status**              |
|:-----------------------------:|
| [![][travis-img]][travis-url] |

The [Discuit.jl][discuit_repo] package contains tools for parameter inference and simulation of discrete state space continuous time (DSSCT) models in Julia. See [Discuit in R][discuit_r_repo] for the corresponding R package.

## Features

- Customisable finite adaptive MCMC algorithm for fast parameter inference.
- Pooley model based proposal (MBP) method for improved mixing.
- Simulation via the Gillespie direct method.
- Automated tools for convergence diagnosis and analysis.

## Installation

The package is not registered and must be added via the package manager Pkg.
From the REPL type `]` to enter the Pkg mode and run:

```
pkg> add https://github.com/mjb3/Discuit.jl
```

## Usage

See the [package documentation][discuit_docs] for further information and examples.

[discuit_repo]: https://github.com/mjb3/Discuit.jl
[discuit_docs]: https://mjb3.github.io/Discuit.jl/latest/
[discuit_r_repo]: https://github.com/mjb3/Discuit
