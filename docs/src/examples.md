# Discuit.jl examples

The following examples provide a flavour of package's core functionality. See the [Discuit.jl manual](@ref) for a description of the data types and functions in Discuit, and [Discuit.jl models](@ref) for a description of the predefined models available in the package.

## Defining a model

`DiscuitModel`s can be created automatically using helper functions or manually by specifying each component. For example the model we are about to create could be generated automatically using `generate_model("SIS", [100,1])`. However constructing it manually is a helpful exercise for getting to know the package. See [Discuit.jl models](@ref) for further details of the `generate_model` function. We start by examining `DiscuitModel` in the package documentation:

```@repl 1
using Discuit;
set_random_seed(1) # hide
?DiscuitModel
```

Next we define a rate function equivalent to that of the basic Kermack-McKendrick `SIS` (and `SIR`) model. The event rates for infection and recovery events respectively are given by:

$r_1 = \theta_1 SI$

$r_2 = \theta_2 I$

The code is correspondingly straightforward:

```@repl 1
function sis_rf(output::Array{Float64, 1}, parameters::Array{Float64, 1}, population::Array{Int64, 1})
    output[1] = parameters[1] * population[1] * population[2]
    output[2] = parameters[2] * population[2]
end
```

Note that the correct signature must be used in the implementation for it to be compatible with the package. Next we define a simple observation function, again with the correct signature:

```@repl 1
obs_fn(population::Array{Int64, 1}) = population
```

The default prior distribution is flat and improper and is equivalent to:

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

model = DiscuitModel("SIS", [100, 1], sis_rf, [-1 1; 1 -1], obs_fn, weak_prior, si_gaussian, 0);
```

## MCMC

The following example is based on that published by Pooley et al. (2015) in the paper that introduces the model based proposal method. The observations data simulated by Pooley can be downloaded [here](https://raw.githubusercontent.com/mjb3/Discuit.jl/master/data/pooley.csv) and saved, e.g. to `path/to/data/`. Next, load the observations data using:

    y = get_observations_from_file("path/to/data/pooley.csv")

Now we can run an MCMC analysis based on the simulated datset:

```@repl 1
y = Observations([20, 40, 60, 80, 100], [0 18; 0 65; 0 70; 0 66; 0 67]); # hide
rs = run_met_hastings_mcmc(model, y, [0.0025, 0.12]);
```

Trace plots are a common visualisation tool in MCMC:

    plot_parameter_trace(rs, 1);

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/jl_traceplot.png" alt="MCMC traceplots" height="220"/>
```

The marginal distribution of parameters can be plotted by calling:

```@repl 1
plot_parameter_marginal(rs, 1)
```

Plotting the data with third party packages such as PyPlot is simple:

    using PyPlot
    x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), parameter]
    plt[:hist](x, 30)

    ylabel("density")

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/marginals.png" alt="MCMC traceplots" height="240"/>
```

A pairwise representation can be produced by calling `plot_parameter_heatmap` (see below for an example).

## Autocorrelation

Autocorrelation can be used to help determine how well the algorithm mixed by using `compute_autocorrelation(rs.mcmc)`. The autocorrelation function for a single Markov chain is implemented in Discuit using the standard formula:

```math
R_l  = \frac{\textrm{E} [(X_i - \bar{X})(X_{i+l} - \bar{X})]}{\sigma^2}
```

for any given lag `l`. The modified formula for multiple chains is given by:

```math
R_{b,l} = \frac{\textrm{E} [ (X_i - \bar{X}_b) ( X_{i + l} - \bar{X}_b ) ]}{\sigma^2_b}
```

$\sigma^2_b = \textrm{E} [(X_i - \bar{X}_b)^2]$

NEED TO ADD image ...

```@repl 1
ac = compute_autocorrelation(rs.mcmc)
p = plot_autocorrelation(ac)
print(p)
```

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/sis-sim.png" alt="SIS simulation" height="180"/>
```

## Convergence diagnostics

### Geweke test of stationarity

The Geweke statistic tests for non-stationarity by comparing the mean and variance for two sections of the Markov chain (see Geweke, 1992; Cowles, 1996). It is given by:

$z = \frac{\bar{\theta}_{i, \alpha} - \bar{\theta}_{i, \beta}}{\sqrt{Var(\theta_{i, \alpha})+Var(\theta_{i, \beta})})}$

```@repl 1
plot_geweke_series(rs)
```

The results of MCMC analyses can also be saved to file for analysis in the companion [R package](https://mjb3.github.io/Discuit/). In Julia, run:

    print_mcmc_results(rs, "path/to/mcmc/data/")

Now, in R, run:

    library(Discuit)
    library(gridExtra)
    rs = LoadMcmcResults("path/to/mcmc/data/")
    tgt = c(0.003, 0.1)
    grid.arrange(PlotGewekeSeries(rs), PlotParameterHeatmap(rs, 1, 2, tgt), nrow = 1)

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/geweke_heatmap.png" alt="MCMC analysis" height="240"/>
```

### Gelman-Rubin convergence diagnostic

The Gelman-Rubin diagnostic is designed to diagnose convergence of two or more Markov chains by comparing within chain variance to between chain variance (Gelman et al, 1992, 2014). The *estimated scale reduction* statistic (sometimes referred to as *potential scale reduction factor*) is calculated for each parameter in the model.

Let ``\bar{\theta}``, $W$ and $B$ be vectors of length $P$ representing the mean of model parameters $\theta$, within chain variance between chain variance respectively for $M$ Markov chains:

$W = \frac{1}{M} \sum_{i = 1}^M \sigma^2_i$

$B = \frac{N}{M - 1} \sum_{i = 1}^M (\hat{\theta}_i - \hat{\theta})^2$

The estimated scale reduction statistic is given by:

$R = \sqrt{\frac{d + 3}{d + 1} \frac{N-1}{N} + (\frac{M+1}{MN} \frac{B}{W})}$

where the first quantity on the RHS adjusts for sampling variance and $d$ is degrees of freedom estimated using the method of moments. For a valid test of convergence the Gelman-Rubin requires two or more Markov chains with over dispersed target values relative to the target distribution. A matrix of such values is therefore required in place of the vector representing the initial values an McMC analysis when calling the function in Discuit, with the $i^{th}$ row vector used to initialise the $i^{th}$ Markov chain.

```@repl 1
rs = run_gelman_diagnostic(model, y, [0.0025 0.08; 0.003 0.12; 0.0035 0.1]);
ac = compute_autocorrelation(rs.mcmc); # hide
```

## Simulation

ADD SIM BLURB `generate_model("LOTKA", [79, 71])`.

```@repl 1
set_random_seed(1); # hide
model = generate_model("LOTKA", [79, 71]); # hide
xi = gillespie_sim(model, [0.5, 0.0025, 0.3]);
```

The simulation can be visualised using `plot_trajectory(xi)`:

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/lotka_sim.png" alt="Lotka Volterra simulation" height="240"/>
```

## Custom MCMC

Some situations...

First we generate a standard `SIR` model and set the `t0_index = 3`.

```@repl 1
set_random_seed(1) # hide
model = generate_model("SIR", [119, 1, 0]);
model.t0_index = 3;
```

Next we define the "medium" prior used by O'Neill and Roberts, with some help from the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package:

```@repl 1
import Pkg; Pkg.add("Distributions"); # hide
using Distributions;
p1 = Gamma(10, 0.0001);
p2 = Gamma(10, 0.01);
function prior_density(parameters::Array{Float64, 1})
    return parameters[3] < 0.0 ? pdf(p1, parameters[1]) * pdf(p2, parameters[2]) * (0.1 * exp(0.1 * parameters[3])) : 0.0
end
model.prior_density = prior_density;
```

The observation model is replaced with one that returns `log(1)` since we will only propose sequences consitent with the observed recoveries and ``\pi(\xi | \theta)`` is evaluated automatically by Discuit):

```@repl 1
observation_model(y::Array{Int, 1}, population::Array{Int, 1}) = 0.0
model.observation_model = observation_model;
```

Next we define an array `t` to contain the recovery times reported by O'Neill and Roberts and a simple `Observations` variable which consists of the maximum event time and an empty two dimensional array:

```@repl 1
t = [0.0, 13.0, 20.0, 22.0, 25.0, 25.0, 25.0, 26.0, 30.0, 35.0, 38.0, 40.0, 40.0, 42.0, 42.0, 47.0, 50.0, 51.0, 55.0, 55.0, 56.0, 57.0, 58.0, 60.0, 60.0, 61.0, 66.0];
y = Observations([67.0], Array{Int64, 2}(undef, 1, 1))
```

We also need to define an initial state using the `generate_custom_x0` function using some parameter values and a vector of event times and corresponding event types, consistent with `t`:

```@repl 1
evt_tm = Float64[];
evt_tp = Int64[];
for i in 1:(length(t) - 1)# infections at arbitrary t > t0
    push!(evt_tm, -4.0)
    push!(evt_tp, 1)
end
for i in eachindex(t)     # recoveries
    push!(evt_tm, t[i])
    push!(evt_tp, 2)
end
x0 = generate_custom_x0(model, y, [0.001, 0.1, -4.0], evt_tm, evt_tp);
```

The final step before we run our analysis is to define the algorithm which will propose changes to augmented data (parameter proposals are automatically configured by Discuit).

```@repl 1
function custom_proposal(model::PrivateDiscuitModel, xi::MarkovState, xf_parameters::ParameterProposal)
    t0 = xf_parameters.value[model.t0_index]
    ## move
    seq_f = deepcopy(xi.trajectory)
    # choose event and define new one
    evt_i = rand(1:length(xi.trajectory.time))
    evt_tm = xi.trajectory.event_type[evt_i] == 1 ? (rand() * (model.obs_data.time[end] - t0)) + t0 : floor(xi.trajectory.time[evt_i]) + rand()
    evt_tp = xi.trajectory.event_type[evt_i]
    # remove old one
    splice!(seq_f.time, evt_i)
    splice!(seq_f.event_type, evt_i)
    # add new one
    if evt_tm > seq_f.time[end]
        push!(seq_f.time, evt_tm)
        push!(seq_f.event_type, evt_tp)
    else
        for i in eachindex(seq_f.time)
            if seq_f.time[i] > evt_tm
                insert!(seq_f.time, i, evt_tm)
                insert!(seq_f.event_type, i, evt_tp)
                break
            end
        end
    end
    # compute ln g(x)
    prop_lk = 1.0
    ## evaluate full likelihood for trajectory proposal and return
    return MarkovState(xi.parameters, seq_f, compute_full_log_like(model, xi.parameters.value, seq_f), prop_lk, 3)
end # end of std proposal function
```

We can now run the MCMC analysis:

```@repl 1
rs = run_custom_mcmc(model, y, custom_proposal, x0, 120000, 20000);
```

Need to add commentary:

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/cmcmc_trace.png" alt="SIR traceplots" height="220"/>
```

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/cmcmc_geweke_hm.png" alt="SIR traceplots" height="220"/>
```


- link to [`set_random_seed(seed::Int64)`](@ref)

## References

```@raw html
<html>
<head>
<script type="text/javascript" src="https://cdn.rawgit.com/pcooksey/bibtex-js/ef59e62c/src/bibtex_js.js"></script>
</head>
<body>
<textarea id="bibtex_input" style="display:none;">
@article{gillespie_exact_1977,
	title = {Exact stochastic simulation of coupled chemical reactions},
	volume = {81},
	issn = {0022-3654, 1541-5740},
	url = {http://pubs.acs.org/doi/abs/10.1021/j100540a008},
	doi = {10.1021/j100540a008},
	language = {en},
	number = {25},
	urldate = {2017-02-18},
	journal = {The Journal of Physical Chemistry},
	author = {Gillespie, Daniel T.},
	month = dec,
	year = {1977},
	pages = {2340--2361}
}

@incollection{geweke_evaluating_1992,
	title = {Evaluating the {Accuracy} of {Sampling}-{Based} {Approaches} to the {Calculation} of {Posterior} {Moments}},
	abstract = {Data augmentation and Gibbs sampling are two closely related, sampling-based approaches to the calculation of posterior moments. The fact that each produces a sample whose constituents are neither independent nor identically distributed complicates the assessment of convergence and numerical accuracy of the approximations to the expected value of functions of interest under the posterior. In this paper methods from spectral analysis are used to evaluate numerical accuracy formally and construct diagnostics for convergence. These methods are illustrated in the normal linear model with informative priors, and in the Tobit-censored regression model.},
	booktitle = {{IN} {BAYESIAN} {STATISTICS}},
	publisher = {University Press},
	author = {Geweke, John},
	year = {1992},
	pages = {169--193}
}

@article{cowles_markov_1996,
	title = {Markov {Chain} {Monte} {Carlo} {Convergence} {Diagnostics}: {A} {Comparative} {Review}},
	volume = {91},
	issn = {01621459},
	shorttitle = {Markov {Chain} {Monte} {Carlo} {Convergence} {Diagnostics}},
	url = {http://www.jstor.org/stable/2291683},
	doi = {10.2307/2291683},
	number = {434},
	urldate = {2018-03-20},
	journal = {Journal of the American Statistical Association},
	author = {Cowles, Mary Kathryn and Carlin, Bradley P.},
	month = jun,
	year = {1996},
	pages = {883}
}

@article{gelman_inference_1992,
	title = {Inference from iterative simulation using multiple sequences},
	journal = {Statistical science},
	author = {Gelman, Andrew and Rubin, Donald B.},
	year = {1992},
	pages = {457--472}
}
@book{gelman_bayesian_2014,
	title = {Bayesian data analysis},
	isbn = {978-1-4398-9820-8 978-1-4398-4096-2},
	language = {English},
	urldate = {2018-03-18},
	author = {Gelman, Andrew and Carlin, John B and Stern, Hal Steven and Dunson, David B and Vehtari, Aki and Rubin, Donald B},
	year = {2014},
	note = {OCLC: 909477393}
}
</textarea>
<div id="bibtex_display"></div>
</body>
</html>
```
