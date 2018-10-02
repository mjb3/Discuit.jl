# Discuit.jl examples

The following examples provide a flavour of Discuit's core functionality. See the [Discuit.jl manual](@ref) for more detailed instructions.

## SIS model

The following example is based on that published by Pooley et al. in 2015 in the paper that introduces the model based proposal method. We could use `generate_model("SIS", [100,1])` to generate the model but constructing it manually is a helpful exercise for getting to know the package. We start by examining `DiscuitModel` in the package documentation:

```@repl 1
using Discuit;
set_random_seed(1) # hide
?DiscuitModel
```

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
model = DiscuitModel("SIS", sis_rf, [-1 1; 1 -1], 0, [100, 1], obs_fn, weak_prior, si_gaussian);
```

### Simulation

Although our main goal is to replicate the analysis of Pooley et al. we can also run a simulation using the `gillespie_sim` function.

```@repl 1
xi = gillespie_sim(model, [0.003, 0.1]);
```

We can also visualise the results using the corresponding R package: rDiscuit. ADD LINK

![SIS simulation](https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/sis-sim.png)

### MCMC

Running an MCMC analysis based on a set of observations data is simple. TBC...

```@repl 1
obs = Observations([20, 40, 60, 80, 100], [0 18; 0 65; 0 70; 0 66; 0 67]);
rs = run_met_hastings_mcmc(model, obs, [0.003, 0.1]);
```

Placeholder for MCMC output.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/traceplots.png" alt="MCMC traceplots" height="240"/>
```

### Diagnostic

#### Geweke

NEED TO ADD geweke definition...

```@repl 1
rs.geweke
```

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/geweke_heatmap.png" alt="MCMC analysis" height="240"/>
```

#### Gelman-Rubin diagnostic

NEED TO ADD gelman definition...

```@repl 1
rs = run_gelman_diagnostic(model, obs, [0.0025 0.08; 0.003 0.12; 0.0035 0.1]);
ac = compute_autocorrelation(rs.mcmc); # hide
```

#### Autocorrelation

Autocorrelation can be used to help determine how well the algorithm mixed by using `compute_autocorrelation(rs.mcmc)`.

NEED TO ADD autocorr definition and image ...

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/Discuit.jl/master/docs/img/sis-sim.png" alt="SIS simulation" height="180"/>
```

## Custom MCMC

Some situations...

First we generate a standard [SIR](@ref) model and set the `t0_index = 3`.

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
