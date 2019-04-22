## resources:
import Discuit
import Random

### pooley SIS model example ###
function pooley()
    Random.seed!(1)
    ## define model
    # rate function
    function sis_rf!(output, parameters::Array{Float64, 1}, population::Array{Int64, 1})
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
    model = Discuit.DiscuitModel("SIS", sis_rf!, [-1 1; 1 -1], [100, 1], obs_fn, prior_density, si_gaussian, 0)

    ## run sim
    xi = Discuit.gillespie_sim(model, [0.003, 0.1])
    println(length(xi.trajectory.time), "\n")
end

function pooley_prebaked()
    Random.seed!(2)

    ## demo 1
    ic = [100, 1]
    model = Discuit.generate_model("SIS", ic);
    theta = [0.003, 0.1]
    x = Discuit.gillespie_sim(model, theta);
    println(Discuit.plot_trajectory(x))

    ## UNCOMMENT ************************
    # ## demo 2
    # y = x.observations
    # initial_theta = [0.005, 0.12]
    # mcmc = run_met_hastings_mcmc(model, x.observations, initial_theta);
    # # NEED TO ADD TAB(MCMC) ***
    #
    # ## demo 3
    # plot_parameter_trace(mcmc, 1)
    # plot_parameter_heatmap(mcmc, 1, 2)


    ## MCMC
    # obs = Observations([20, 40, 60, 80, 100], [0 18; 0 65; 0 70; 0 66; 0 67]);
    # obs = Discuit.get_observations("./data/pooley.csv")
    obs = x.observations
    rs = Discuit.run_met_hastings_mcmc(model, obs, [0.003, 0.1]);
    Discuit.tabulate_mcmc_results(rs, true)
    ac = Discuit.compute_autocorrelation(rs)
    # # print
    # print_mcmc_results(rs, "./out/doc/mcmc_example/")
    #
    # ## Diagnostics
    # # geweke
    # println(" geweke statistics: ", rs.geweke[2][1,:], "\n")
    # gelman
    rs = Discuit.run_gelman_diagnostic(model, obs, [0.002 0.08; 0.0028 0.12; 0.0035 0.1])
    Discuit.tabulate_gelman_results(rs, true)
    # print_gelman_results(rs, "./out/gelman_example/")
    # # # autocorrelation
    # ac = compute_autocorrelation(rs.mcmc)
    # print_autocorrelation(ac, string("./out/doc/acp_mbp.csv"))
    #
    # standard proposals (for comparison)
    rs = Discuit.run_gelman_diagnostic(model, obs, [0.0025 0.08; 0.003 0.12; 0.0035 0.1], 200000, 20000, false);
    Discuit.tabulate_gelman_results(rs, true)
    println(Discuit.plot_parameter_trace(rs.mcmc, 1))
    # ac = compute_autocorrelation(rs.mcmc)
    # print_autocorrelation(ac, string("./out/doc/acp_std.csv"))
end

## custom roberts
import Distributions
function custom_bobs()
    # set_random_seed(1)
    ## define model
    model = Discuit.generate_model("SIR", [119, 1, 0]);
    model.t0_index = 3
    # add "medium" prior
    p1 = Distributions.Gamma(10, 0.0001)
    p2 = Distributions.Gamma(10, 0.01)
    function prior_density(parameters::Array{Float64, 1})
        return parameters[3] < 0.0 ? Distributions.pdf(p1, parameters[1]) * Distributions.pdf(p2, parameters[2]) * (0.1 * exp(0.1 * parameters[3])) : 0.0
    end
    # 'weak' prior
    function prior_density(parameters::Array{Float64, 1})
        return parameters[3] < 0.0 ? Distributions.pdf(p1, parameters[1]) * Distributions.pdf(p2, parameters[2]) * (0.1 * exp(0.1 * parameters[3])) : 0.0
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
    x0 = Discuit.generate_custom_x0(model, y, [0.001, 0.1, -4.0], evt_tm, evt_tp);
    println("x0 log like: ", x0.log_like)

    ## custom proposal algorithm
    # - IS THERE A MORE EFFICIENT WAY TO DO THIS? I.E. ROTATE using circshift or something?
    function custom_proposal(model::Discuit.PrivateDiscuitModel, xi::Discuit.MarkovState, xf_parameters::Discuit.ParameterProposal)
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
        return Discuit.MarkovState(xi.parameters, seq_f, Discuit.compute_full_log_like(model, xi.parameters.value, seq_f), prop_lk, 3)
    end # end of std proposal function

    ## run MCMC
    rs = run_custom_mcmc(model, y, custom_proposal, x0, 120000, 20000)
    Discuit.tabulate_mcmc_results(rs, true)
    # print_mcmc_results(rs, "./out/doc/custom_mcmc_example/")
end

# plot_parameter_trace
# pooley()
pooley_prebaked()
# custom_bobs()
