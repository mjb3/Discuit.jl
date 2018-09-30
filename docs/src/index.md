# Discuit.jl documentation

*Fast parameter inference for discrete state space continuous time (DSSCT) models in Julia.*
<!-- Discuit: simulation and parameter inference for discrete state space continuous time (DSSCT) models. -->

!!! note
    Please note that this package is still in development.


## Package Features

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

## Usage

The [Discuit.jl examples](@ref) section provides enough code to get up and running. A more detailed guide to the pre defined models is available in the [Discuit.jl models](@ref) section. Further information regarding the packages other functionality can be found in the [Discuit.jl manual](@ref).

## test

The autocorrelation function for a single Markov chain is implemented in Discuit using the standard formula:

```math
R_l  = \frac{\ev[(X_i - \bar{X})(X_{i+l} - \bar{X})]}{\sigma^2}
```

for any given lag `l`. The modified formula for multiple chains is given by:

```math
R^{\prime}_l  = \frac{\ev[(X_i - \bar{X}_b)(X_{i+l} - \bar{X}_b)]}{\sigma^2_b}\\
```

```math
\sigma^2_b = \ev[(X_i - \bar{X}_b)^2]
```
