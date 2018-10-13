using Discuit
# import PyPlot
import UnicodePlots

function traj_example()
    set_random_seed(1)
    model = generate_model("LOTKA", [79, 71]);

    ## run sim
    x = gillespie_sim(model, [0.5, 0.0025, 0.3]);
    # plot_trajectory(xi)
    p = UnicodePlots.lineplot(x.trajectory.time, x.population[:,1], title = string(x.model_name, " simulation"), name = string(x.model_name[1]))
    for i in 2:size(x.population, 2)
        UnicodePlots.lineplot!(p, x.trajectory.time, x.population[:,i], name = string(x.model_name[i]))
    end
    UnicodePlots.xlabel!(p, "time")
    UnicodePlots.ylabel!(p, "population")
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
    rs = run_met_hastings_mcmc(model, y, [0.003, 0.1])
    ## traceplots
    # PyPlot.subplot(1, 2, 1)
    # plot_parameter_trace(rs, 1)


    x = mcmc.adapt_period:size(mcmc.samples, 1)
    PyPlot.plot(x, mcmc.samples[mcmc.adapt_period:size(mcmc.samples, 1), parameter])
    PyPlot.xlabel("sample")
    PyPlot.ylabel(string("\$\\theta_", parameter, "\$"))
    # PyPlot.subplot(1, 2, 2)
    # plot_parameter_trace(rs, 2)
    # PyPlot.show()
    ## marginal
    # PyPlot.subplot(1, 2, 1)
    # plot_parameter_marginal(rs, 1)
    # PyPlot.subplot(1, 2, 2)
    # plot_parameter_marginal(rs, 2)
    # PyPlot.show()
    ## heatmap
    # plot_parameter_heatmap(rs, 1, 2)
    ## geweke plot
    x = rs.geweke[1]
    print("x: ", x, "\n")
    PyPlot.scatter(x, rs.geweke[2][:,1])
    PyPlot.scatter(x, rs.geweke[2][:,1])
end

traj_example()
# mcmc_example()
