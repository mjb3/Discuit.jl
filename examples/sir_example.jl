## resources:
import Discuit
import Random

## example
function sir_example()
    Random.seed!(1)
    # generate model
    initial_condition = [50, 1, 0]
    model = Discuit.generate_model("SIR", initial_condition);

    ## sim
    theta = [0.005, 0.1]
    x = Discuit.gillespie_sim(model, theta, 100.0, 10);
    # Discuit.print_trajectory(model, x, "./out/sir_sim.csv")
    println(Discuit.plot_trajectory(x))

    obs = x.observations
    theta_init = [0.001 0.1; 0.002 0.1; 0.003 0.1; 0.001 0.2; 0.002 0.2; 0.003 0.2]
    # Discuit.print_observations(obs, "./out/sir_obs.csv")

    Random.seed!(6)
    ## MCMC (STD)
    # rs = Discuit.run_met_hastings_mcmc(model, obs, [0.001, 0.1], 300000, 50000, false)
    # Discuit.tabulate_mcmc_results(rs, true)
    # println(Discuit.plot_parameter_trace(rs, 1))
    # println(Discuit.plot_parameter_heatmap(rs, 1, 2))
    ## convergence
    gmn = Discuit.run_gelman_diagnostic(model, obs, theta_init, 100000, 20000, false)
    Discuit.tabulate_gelman_results(gmn, true)

    ## MCMC (MBP)
    # rs = Discuit.run_met_hastings_mcmc(model, obs, [0.001, 0.1], 100000, 20000)
    # Discuit.tabulate_mcmc_results(rs, true)
    # println(Discuit.plot_parameter_trace(rs, 1))
    # println(Discuit.plot_parameter_heatmap(rs, 1, 2))
    gmn = Discuit.run_gelman_diagnostic(model, obs, theta_init)
    Discuit.tabulate_gelman_results(gmn, true)

end
sir_example()
