using Discuit
# import PyPlot
import UnicodePlots

function traj_example()
    set_random_seed(1)
    model = generate_model("LOTKA", [79, 71]);

    ## run sim
    x = gillespie_sim(model, [0.5, 0.0025, 0.3]);
    p = plot_trajectory(x)
    # p = UnicodePlots.lineplot(x.trajectory.time, x.population[:,1], title = string(x.model_name, " simulation"), name = string(x.model_name[1]))
    # for i in 2:size(x.population, 2)
    #     UnicodePlots.lineplot!(p, x.trajectory.time, x.population[:,i], name = string(x.model_name[i]))
    # end
    # UnicodePlots.xlabel!(p, "time")
    # UnicodePlots.ylabel!(p, "population")
    print(p)
end

import HTTP
import DelimitedFiles

function mcmc_example()
    set_random_seed(1)
    model = generate_model("SIS", [100, 1]);
    ## download pooley obs
    res = HTTP.get("https://raw.githubusercontent.com/mjb3/Discuit.jl/master/data/pooley.csv")
    df = DelimitedFiles.readdlm(res.body, ','; header=true)
    y = get_observations_from_array(df[1])

    ## mcmc
    mcmc = run_met_hastings_mcmc(model, y, [0.003, 0.1])
    ## traceplots
    # PyPlot.subplot(1, 2, 1)

    # p = plot_parameter_trace(mcmc, 1)
    # print(p)

    # x = mcmc.adapt_period:size(mcmc.samples, 1)
    # p = UnicodePlots.lineplot(x, mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1)], title = string("F₁", " traceplot."))
    # UnicodePlots.xlabel!(p, "sample")
    # UnicodePlots.ylabel!(p, string("θ", Char(8320 + 1)))
    # PyPlot.subplot(1, 2, 2)
    # plot_parameter_trace(rs, 2)
    # PyPlot.show()

    ## marginal
    # PyPlot.subplot(1, 2, 1)
    # plot_parameter_marginal(rs, 1)

    # parameter = 1
    # x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), parameter]
    # p = UnicodePlots.histogram(x, bins = 20)
    # UnicodePlots.ylabel!(p, string("θ", Char(8320 + parameter)))
    # UnicodePlots.xlabel!(p, "density")
    # print(p)

    # PyPlot.subplot(1, 2, 2)
    # plot_parameter_marginal(rs, 2)
    # PyPlot.show()
    ## heatmap
    # plot_parameter_heatmap(rs, 1, 2)
    x = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), 1]
    y = mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), 2]
    p = UnicodePlots.densityplot(x, y, color = :red)
    print(p)
    ## geweke plot
    # x = rs.geweke[1]
    # print("x: ", x, "\n")
    # PyPlot.scatter(x, rs.geweke[2][:,1])
    # PyPlot.scatter(x, rs.geweke[2][:,1])
end

traj_example()
mcmc_example()
