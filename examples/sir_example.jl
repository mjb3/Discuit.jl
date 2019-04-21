## resources:
import Discuit
import Random

## example
function sir_example()
    ## model
    ic = [100, 1, 0]
    model = Discuit.generate_model("SIR", ic);
    ## sim
    theta = [0.002, 0.1]
    x = Discuit.gillespie_sim(model, theta);
    println(Discuit.plot_trajectory(x))

end
sir_example()
