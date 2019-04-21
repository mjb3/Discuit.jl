## resources:
import Discuit
import Random

## example
function sir_example()
    Random.seed!(1)
    ## model
    ic = [100, 1, 0]
    model = Discuit.generate_model("SIR", ic);

    ## sim
    theta = [0.002, 0.1]
    x = Discuit.gillespie_sim(model, theta);
    println(Discuit.plot_trajectory(x))

    ## MCMC (MBP)
    obs = x.observations
    rs = Discuit.run_met_hastings_mcmc(model, obs, [0.001, 0.1])
    Discuit.tabulate_mcmc_results(rs, true, 100000, 20000)
    println(Discuit.plot_parameter_trace(rs, 1))
    println(Discuit.plot_parameter_heatmap(rs, 1, 2))

    ## MCMC (STD)
    obs = x.observations
    rs = Discuit.run_met_hastings_mcmc(model, obs, [0.001, 0.1], 300000, 50000, false)
    Discuit.tabulate_mcmc_results(rs, true)
    println(Discuit.plot_parameter_trace(rs, 1))
    println(Discuit.plot_parameter_heatmap(rs, 1, 2))

end
sir_example()
