using Discuit
# import PyPlot
import UnicodePlots

function traj_example()
    set_random_seed(1)
    model = generate_model("LOTKA", [79, 71]);

    ## run sim
    x = gillespie_sim(model, [0.5, 0.0025, 0.3]);
    p = plot_trajectory(x)
    println(p)
end

import HTTP
import DelimitedFiles

function mcmc_example()
    set_random_seed(1)
    model = generate_model("SIS", [100, 1]);
    ## download pooley obs
    res = HTTP.get("https://raw.githubusercontent.com/mjb3/Discuit.jl/master/data/pooley.csv")
    df = DelimitedFiles.readdlm(res.body, ','; header=true)
    y = get_observations(df[1])

    ## mcmc
    mcmc = run_met_hastings_mcmc(model, y, [0.003, 0.1])
    ## traceplots
    # p = plot_parameter_trace(mcmc, 1)
    # print(p)

    ## marginal
    # p = plot_parameter_marginal(mcmc, 1)
    # print(p)

    ## heatmap
    # p = plot_parameter_heatmap(mcmc, 1, 2)
    # print(p)

    ## geweke plot
    p = plot_geweke_series(mcmc)
    print(p)

    # x = mcmc.geweke[1]
    # # print("x: ", x, "\n")
    # p = UnicodePlots.scatterplot(x, mcmc.geweke[2][:,1])
    # UnicodePlots.scatterplot!(p, x, mcmc.geweke[2][:,2])
    # UnicodePlots.lineplot!(p, -2.0, 0.0, color = :yellow)
    # UnicodePlots.lineplot!(p, 2.0, 0.0, color = :yellow)

    ## autocorrelation
    # ac = compute_autocorrelation(mcmc)
    # p = plot_autocorrelation(ac)
    # print(p)
end

# traj_example()
mcmc_example()
