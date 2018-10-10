using Discuit
# import PyPlot

function traj_example()
    set_random_seed(1)
    model = generate_model("LOTKA", [79, 71]);

    ## run sim
    xi = gillespie_sim(model, [0.5, 0.0025, 0.3]);
    plot_trajectory(xi)
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
    ## traceplot
    plot_parameter_trace(rs, 1)
    # plot_parameter_trace(rs, 2)
    ## marginal
    # plot_parameter_marginal(rs, 1)
    # plot_parameter_marginal(rs, 2)
    ## heatmap
    # plot_parameter_heatmap(rs, 1, 2)
end

traj_example()
# mcmc_example()
