# Discuit.jl documentation

*Fast parameter inference for discrete state space continuous time (DSSCT) models in Julia.*
<!-- Discuit: simulation and parameter inference for discrete state space continuous time (DSSCT) models. -->

!!! note
    Please note that this package is still in development.

## Package Features

```@docs
Discuit
```

- User defined DSSCT models.
- Pre programmed with many well known epidemiological models.
- Exact simulation using Gillespie's algorithm.
- Data augmented Markov chain Monte Carlo (MCMC).
- Automated autocorrelation; Geweke and Gelman-Rubin diagnostics.
- Developed for Julia `1.0`.

## Contents

```@contents
```

## Installation

The package can be installed by typing `]` in the REPL to enter the Pkg mode and running:

```
pkg> add https://github.com/mjb3/Discuit.jl
```

## Usage

The [Discuit.jl examples](@ref) section provides enough code to get up and running with further information available in the [Discuit.jl manual](@ref).

## Index

```@index
```
