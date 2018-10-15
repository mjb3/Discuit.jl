## import the package
using Discuit

### pooley ###
function pooley()
    set_random_seed(1)
    ## define model
    # rate function
    function sis_rf(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
        output[1] = parameters[1] * population[1] * population[2]
        output[2] = parameters[2] * population[2]
    end
    # define obs function (no error)
    obs_fn(population::Array{Int64, 1}) = population
    # prior
    function prior_density(parameters::Array{Float64, 1})
        parameters[1] > 0.0 || return 0.0
        parameters[2] > 0.0 || return 0.0
        return 1.0
    end
    # obs model
    function si_gaussian(y::Array{Int64, 1}, population::Array{Int64, 1})
        obs_err = 2
        tmp1 = log(1 / (sqrt(2 * pi) * obs_err))
        tmp2 = 2 * obs_err * obs_err
        obs_diff = y[2] - population[2]
        return tmp1 - ((obs_diff * obs_diff) / tmp2)
    end
    # define model
    model = DiscuitModel("SIS", sis_rf, [-1 1; 1 -1], [100, 1], obs_fn, prior_density, si_gaussian, 0)

    ## run sim
    xi = gillespie_sim(model, [0.003, 0.1])
    print(length(xi.trajectory.time), "\n")
end

function pooley_prebaked()
    set_random_seed(1)
    model = generate_model("SIS", [100, 1]);

    ## run sim
    xi = gillespie_sim(model, [0.003,0.1]);

    ## MCMC
    obs = Observations([20, 40, 60, 80, 100], [0 18; 0 65; 0 70; 0 66; 0 67]);
    # print_observations(obs, "./docs/data/pooley.csv")
    # rs = run_met_hastings_mcmc(model, obs, [0.003, 0.1]);
    # # print
    # print_mcmc_results(rs, "./out/doc/mcmc_example/")
    #
    # ## Diagnostics
    # # geweke
    # print(" geweke statistics: ", rs.geweke[2][1,:], "\n")
    # # gelman
    rs = run_gelman_diagnostic(model, obs, [0.0025 0.08; 0.003 0.12; 0.0035 0.1]);
    print_gelman_results(rs, "./out/gelman_example/")
    # # # autocorrelation
    # ac = compute_autocorrelation(rs.mcmc)
    # print_autocorrelation(ac, string("./out/doc/acp_mbp.csv"))
    #
    # # standard proposals (for comparison)
    # rs = run_gelman_diagnostic(model, obs, [0.0025 0.08; 0.003 0.12; 0.0035 0.1], 80000, 30000, false);
    # ac = compute_autocorrelation(rs.mcmc)
    # print_autocorrelation(ac, string("./out/doc/acp_std.csv"))
end

## custom roberts
using Distributions
function custom_bobs()
    set_random_seed(1)
    ## define model
    model = generate_model("SIR", [119, 1, 0]);
    model.t0_index = 3
    # add "medium" prior
    p1 = Gamma(10, 0.0001)
    p2 = Gamma(10, 0.01)
    function prior_density(parameters::Array{Float64, 1})
        return parameters[3] < 0.0 ? pdf(p1, parameters[1]) * pdf(p2, parameters[2]) * (0.1 * exp(0.1 * parameters[3])) : 0.0
    end
    # 'weak' prior
    function prior_density(parameters::Array{Float64, 1})
        return parameters[3] < 0.0 ? pdf(p1, parameters[1]) * pdf(p2, parameters[2]) * (0.1 * exp(0.1 * parameters[3])) : 0.0
    end
    model.prior_density = prior_density
    # dummy observation model
    observation_model(y::Array{Int, 1}, population::Array{Int, 1}) = 0.0
    model.observation_model = observation_model

    ## initial trajectory
    # removal times
    t = [0.0, 13.0, 20.0, 22.0, 25.0, 25.0, 25.0, 26.0, 30.0, 35.0, 38.0, 40.0, 40.0, 42.0, 42.0, 47.0, 50.0, 51.0, 55.0, 55.0, 56.0, 57.0, 58.0, 60.0, 60.0, 61.0, 66.0];
    y = Observations([67.0], Array{Int64, 2}(undef, 1, 1));
    # initial sequence
    # n::Int64 = (2 * length(t)) - 1;
    evt_tm = Float64[];
    evt_tp = Int64[];
    # infections ar arbitrary t (must be > t0)
    for i in 1:(length(t) - 1)
        push!(evt_tm, -4.0)
        push!(evt_tp, 1)
    end
    # recoveries
    for i in eachindex(t)
        push!(evt_tm, t[i])
        push!(evt_tp, 2)
    end
    x0 = generate_custom_x0(model, y, [0.001, 0.1, -4.0], evt_tm, evt_tp);
    print("x0 log like: ", x0.log_like, "\n")

    ## custom proposal algorithm
    # - IS THERE A MORE EFFICIENT WAY TO DO THIS? I.E. ROTATE using circshift or something?
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

    ## run MCMC
    rs = run_custom_mcmc(model, y, custom_proposal, x0, 120000, 20000)
    print_mcmc_results(rs, "./out/doc/custom_mcmc_example/")
end


# pooley()
pooley_prebaked()
# custom_bobs()
