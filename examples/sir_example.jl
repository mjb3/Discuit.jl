## resources:
import Discuit
import Random

## example
function sir_example()
    Random.seed!(1)
    ## model
    ic = [50, 1, 0]
    model = Discuit.generate_model("SIR", ic);

    ## sim
    theta = [0.005, 0.1]
    x = Discuit.gillespie_sim(model, theta, 100.0, 10);
    # Discuit.print_trajectory(model, x, "./out/sir_sim.csv")
    println(Discuit.plot_trajectory(x))

    obs = x.observations
    # Discuit.print_observations(obs, "./out/sir_obs.csv")

    Random.seed!(1)
    ## MCMC (STD)
    # rs = Discuit.run_met_hastings_mcmc(model, obs, [0.001, 0.1], 300000, 50000, false)
    # Discuit.tabulate_mcmc_results(rs, true)
    # println(Discuit.plot_parameter_trace(rs, 1))
    # println(Discuit.plot_parameter_heatmap(rs, 1, 2))
    ## convergence
    gmn = Discuit.run_gelman_diagnostic(model, obs, [0.001 0.1; 0.002 0.1; 0.003 0.1], 100000, 20000, false)
    Discuit.tabulate_gelman_results(gmn, true)

    ## MCMC (MBP)
    # rs = Discuit.run_met_hastings_mcmc(model, obs, [0.001, 0.1], 100000, 20000)
    # Discuit.tabulate_mcmc_results(rs, true)
    # println(Discuit.plot_parameter_trace(rs, 1))
    # println(Discuit.plot_parameter_heatmap(rs, 1, 2))
    gmn = Discuit.run_gelman_diagnostic(model, obs, [0.001 0.1; 0.002 0.1; 0.003 0.1])
    Discuit.tabulate_gelman_results(gmn, true)

end
sir_example()
