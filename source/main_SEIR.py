from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint

from source.model import number_of_E_compartments, number_of_I_compartments, EpidemicModel
from source.utils import print_endpoint, plot_cumulative, plot_all, plot_final_sizes, \
    plot_and_save_all, plot_and_save_compare

plt.style.use('ggplot')
chn_population = 1436980159


def main():
    # Booleans for turning on/off cases
    bool_solve_orig_seir = True
    bool_solve_controlled_seir = True

    # Solve original SEIR
    init_for_control = []
    if bool_solve_orig_seir:
        init_for_control = solve_orig_seir()

    # Solve controlled SEIR
    if bool_solve_controlled_seir:
        solve_controlled_seir(init_for_control)


def solve_orig_seir():
    # Set rates
    incubation_period = 5
    infectious_period = 10
    r_0 = 2.6
    beta = r_0 / (chn_population * infectious_period)
    alpha = 1 / (incubation_period / number_of_E_compartments)
    gamma = 1 / (infectious_period / number_of_I_compartments)
    params = np.array((beta, alpha, gamma))
    # Set initial values
    init_values = {"s0": [chn_population - 41],
                   "e0": [41, 0],
                   "i0": [0, 0, 0],
                   "r0": [0],
                   "c0": [0]
                   }
    # Instantiate epidemic model
    model = EpidemicModel(init_values)
    x0 = np.array(model.get_initial_values())
    # Solve epidemic model
    t = np.linspace(0, 32, 100000)
    solution = solve_model(t, x0, params, model)
    cumulative = solution[:, -1]
    # Plot cumulative compartment
    plot_cumulative(cumulative, t)
    # Plot all compartments
    plot_all(solution, t)
    # Get state value corresponding to specified cumulative state
    init_for_control = print_endpoint(solution, cumulative, t)
    return init_for_control


def solve_controlled_seir(init_for_control):
    # Hubei population
    hubei_population = 59170000
    # Hopkins: 23.01.2020. total number of infected = 865, infected cases NOT from Hubei province = 316
    outside_hubei_ratio = 316 / 865
    init = outside_hubei_ratio * init_for_control
    s0 = chn_population - (hubei_population - np.sum((1 - outside_hubei_ratio) * init_for_control[1:-1]))
    # Set time intervals
    incubation_period = 5
    infectious_period = 10
    # List of R0s
    r_0_list = [1.95, 2.15, 2.13, 2.2, 2.6, 2.9]
    # Reported cases after 23. January (Hopkins)
    cases_after_23jan = {"0": 316, "0.5": 367, "1": 591,
                         "1.5": 638, "2": 1004, "2.5": 1314, "3": 1402, "3.5": 1815, "4": 2503, "5": 2610}
    # Initial list for best stars determined per R0
    best_t_stars = []
    # Initial list for final sizes determined per R0 and t*
    final_sizes = []
    for r_0 in r_0_list:
        # Set rates
        beta = r_0 / ((chn_population - hubei_population) * infectious_period)
        alpha = 1 / (incubation_period / number_of_E_compartments)
        gamma = 1 / (infectious_period / number_of_I_compartments)
        params = np.array((beta, alpha, gamma))
        # Set initial values
        init_values = {"s0": [s0],
                       "e0": [init[1], init[2]],
                       "i0": [init[3], init[4], init[5]],
                       "r0": [init[6]],
                       "c0": [init[7]]
                       }
        # Initial list for solutions for all t_star
        sols = []
        # Initial list for final sizes for all t_star
        fs = []
        for t_star in np.arange(10, 51, 1):
            # Instantiate epidemic model
            model = EpidemicModel(init_values, t_star=t_star, r_0=r_0)
            x0 = np.array(model.get_initial_values())
            # Solve epidemic model
            t = np.linspace(0, 10 * t_star, 100000)
            solution = solve_model(t, x0, params, model)
            # Save solution for current t_star
            sols.append(solution)
            # Save final size for current t_star
            fs.append(solution[-1, -1])

            # Transform reported cases data to lists
            t_data = np.array(list(cases_after_23jan.keys())).astype(float)
            cases_data = np.array(list(cases_after_23jan.values())).astype(float)
            # Plot solutions and save figures
            plot_and_save_all(t, solution, r_0, t_star)
            # Plot solutions and data and save figures
            plot_and_save_compare(cases_data, r_0, solution, t, t_data, t_star)

        # Find t_star, which corresponds to solution fitting best to data
        best_t_star = find_closest_solution(sols, cases_after_23jan)
        # Store best t_star for current R0
        best_t_stars.append(best_t_star)
        # Store final sizes for current R0
        final_sizes.append(fs)
    # Print best t_stars
    for (r_0, t_star) in zip(r_0_list, best_t_stars):
        print("Best t* for R0=", r_0, "is:", t_star)
    # Plot and save t_star vs. final size figures for all R0s
    for (r_0, final_size) in zip(r_0_list, final_sizes):
        plot_final_sizes(r_0, final_size)


def find_closest_solution(solution, cases: dict):
    """
    Find best t_star value via squared error between data and solution
    """
    # Smallest t_star
    t_star = 10
    # Init values for minimum search
    t_min = -1
    sse_min = 10e9
    for sol in solution:
        sum_diff_square = 0
        for t in cases.keys():
            index = math.ceil(100000 * float(t) / (10 * t_star))
            sum_diff_square += (sol[index, -1] - cases[t]) ** 2
        if sum_diff_square < sse_min:
            sse_min = sum_diff_square
            t_min = t_star
        t_star += 1
    return t_min


def solve_model(t, x0, params, model):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(model.get_model, x0, t, args=(params,))
    return x


if __name__ == "__main__":
    main()
