from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from source.model import EpidemicModel
from source.utils import plot_final_sizes, plot_and_save_all

plt.style.use('ggplot')


def main():
    # Booleans for turning on/off cases
    bool_solve_controlled_seir = True

    # Solve controlled SEIR
    if bool_solve_controlled_seir:
        solve_controlled_seir()


def solve_controlled_seir():
    # Hubei & China population
    hubei_population = 58160000
    chn_population = 1437000000
    population = chn_population - hubei_population
    # Set time intervals
    incubation_period = 4.8
    # List of R0s
    r_0_list = [2.1, 2.6, 3.1]
    # List of infectious periods associated to R0s
    infectious_period_list = [2.8, 4, 5.1]
    # Initial list for final sizes determined per R0 and t*
    final_sizes = []
    # List of t*-s
    t_stars = np.arange(20, 61, 1)
    for (r_0, infectious_period) in zip(r_0_list, infectious_period_list):
        # Set rates
        beta = r_0 / (population * infectious_period)
        alpha = 1 / incubation_period
        gamma = 1 / infectious_period
        params = np.array((beta, alpha, gamma))
        # Set initial values
        init_values = {"e0": [6000, 7500],
                       "i0": [1750, 1560, 1413],
                       "r0": [0],
                       "c0": [4723]
                       }
        not_susceptible = sum([item for sublist in init_values.values() for item in sublist])
        init_values.update({"s0": [population - not_susceptible]})
        # Initial list for solutions for all t_star
        sols = []
        # Initial list for final sizes for all t_star
        fs = []
        for t_star in t_stars:
            # Instantiate epidemic model
            model = EpidemicModel(init_values, t_star=t_star, r_0=r_0)
            x0 = np.array(model.get_initial_values())
            # Solve epidemic model
            t = np.linspace(0, 150, 100000)
            solution = solve_model(t, x0, params, model)
            # Save solution for current t_star
            sols.append(solution)
            # Save final size for current t_star
            fs.append(solution[-1, -1])

            # Plot solutions and save figures
            plot_and_save_all(t, solution, r_0, t_star)

        # Store final sizes for current R0
        final_sizes.append(fs)
    # Plot and save t_star vs. final size figures for all R0s
    for (r_0, final_size) in zip(r_0_list, final_sizes):
        plot_final_sizes(r_0, t_stars, final_size)


def solve_model(t, x0, params, model):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(model.get_model, x0, t, args=(params,))
    return x


if __name__ == "__main__":
    main()
