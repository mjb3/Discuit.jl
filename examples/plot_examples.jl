using Discuit
using PyPlot

function traj_example()
    set_random_seed(1)
    model = generate_model("LOTKA", [79, 71]);

    ## run sim
    xi = gillespie_sim(model, [0.003,0.1]);

    #
    plot(xi.trajectory.time, xi.population)
end

function mcmc_example()
    set_random_seed(1)
    model = generate_model("SIS", [100, 1]);

    ## tbc
end

traj_example()
