using Discuit
import PyPlot

function traj_example()
    set_random_seed(1)
    model = generate_model("LOTKA", [79, 71]);

    ## run sim
    xi = gillespie_sim(model, [0.5, 0.0025, 0.3]);

    #
    PyPlot.plot(xi.trajectory.time, xi.population)
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
    x = rs.adapt_period:size(rs.samples, 1)
    y = rs.samples[rs.adapt_period:size(rs.samples, 1), 2]
    PyPlot.plot(x, y)
    ## marginal
end

# traj_example()
mcmc_example()
