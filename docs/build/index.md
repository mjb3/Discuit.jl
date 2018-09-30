
<a id='Discuit.jl-documentation-1'></a>

# Discuit.jl documentation


*Fast parameter inference for discrete state space continuous time (DSSCT) models in Julia.* <!– Discuit: simulation and parameter inference for discrete state space continuous time (DSSCT) models. –>


!!! note
    Please note that this package is still in development.



The seminal work [see @pizza2000identification]


<a id='Package-Features-1'></a>

## Package Features

<a id='Discuit' href='#Discuit'>#</a>
**`Discuit`** &mdash; *Module*.



**module Discuit**

Discuit is a package for:

  * User defined DSSCT models.
  * Pre programmed with many well known epidemiological models.
  * Exact simulation using Gillespie's algorithm.
  * Data augmented Markov chain Monte Carlo (MCMC).
  * Automated autocorrelation; Geweke and Gelman-Rubin diagnostics.
  * Developed for Julia `1.0`.
  * Author: Martin Burke (martin.burke@bioss.ac.uk)
  * Date: 2018-08-22


<a target='_blank' href='https://github.com/mjb3/Discuit.jl/blob/05edd537241a52c08fc338fccfabed0c6ba0302f/src/Discuit.jl#L1-L15' class='documenter-source'>source</a><br>


<a id='Contents-1'></a>

## Contents

- [Discuit.jl models](models.md#Discuit.jl-models-1)
    - [Model builder](models.md#Model-builder-1)
    - [Classic Kermack-McKendrick models](models.md#Classic-Kermack-McKendrick-models-1)
    - [Latent Kermack-McKendrick models](models.md#Latent-Kermack-McKendrick-models-1)
    - [Miscellaneous](models.md#Miscellaneous-1)
- [Discuit.jl examples](examples.md#Discuit.jl-examples-1)
    - [SIS model](examples.md#SIS-model-1)
    - [Custom MCMC](examples.md#Custom-MCMC-1)
- [Discuit.jl manual](manual.md#Discuit.jl-manual-1)
    - [Contents](manual.md#Contents-1)
    - [Custom structs](manual.md#Custom-structs-1)
    - [Functions](manual.md#Functions-1)
    - [Index](manual.md#Index-1)
- [Discuit.jl documentation](index.md#Discuit.jl-documentation-1)
    - [Package Features](index.md#Package-Features-1)
    - [Contents](index.md#Contents-1)
    - [Installation](index.md#Installation-1)
    - [Usage](index.md#Usage-1)
    - [Index](index.md#Index-1)


<a id='Installation-1'></a>

## Installation


The package can be installed by typing `]` in the REPL to enter the Pkg mode and running:


```
pkg> add https://github.com/mjb3/Discuit.jl
```


<a id='Usage-1'></a>

## Usage


The [Discuit.jl examples](examples.md#Discuit.jl-examples-1) section provides enough code to get up and running with further information available in the [Discuit.jl manual](manual.md#Discuit.jl-manual-1).


<a id='Index-1'></a>

## Index

- [`Discuit`](index.md#Discuit)
- [`Discuit.DiscuitModel`](examples.md#Discuit.DiscuitModel)
- [`Discuit.GelmanResults`](manual.md#Discuit.GelmanResults)
- [`Discuit.MCMCResults`](manual.md#Discuit.MCMCResults)
- [`Discuit.Observations`](manual.md#Discuit.Observations)
- [`Discuit.SimResults`](manual.md#Discuit.SimResults)
- [`Discuit.generate_gaussian_obs_model`](manual.md#Discuit.generate_gaussian_obs_model)
- [`Discuit.generate_generic_obs_function`](manual.md#Discuit.generate_generic_obs_function-Tuple{})
- [`Discuit.generate_model`](manual.md#Discuit.generate_model)
- [`Discuit.generate_weak_prior`](manual.md#Discuit.generate_weak_prior-Tuple{Int64})
- [`Discuit.gillespie_sim`](manual.md#Discuit.gillespie_sim)
- [`Discuit.print_gelman_results`](manual.md#Discuit.print_gelman_results-Tuple{GelmanResults,String})
- [`Discuit.print_mcmc_results`](manual.md#Discuit.print_mcmc_results-Tuple{MCMCResults,String})
- [`Discuit.print_observations`](manual.md#Discuit.print_observations-Tuple{Observations,String})
- [`Discuit.print_trajectory`](manual.md#Discuit.print_trajectory-Tuple{DiscuitModel,SimResults,String})
- [`Discuit.read_obs_data_from_file`](manual.md#Discuit.read_obs_data_from_file-Tuple{String})
- [`Discuit.run_custom_mcmc`](manual.md#Discuit.run_custom_mcmc)
- [`Discuit.run_gelman_diagnostic`](manual.md#Discuit.run_gelman_diagnostic)
- [`Discuit.run_met_hastings_mcmc`](manual.md#Discuit.run_met_hastings_mcmc)
- [`Discuit.set_random_seed`](manual.md#Discuit.set_random_seed-Tuple{Int64})

