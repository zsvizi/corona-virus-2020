from __future__ import division
import numpy as np
from scipy.integrate import odeint, solve_ivp

from source.model import EpidemicModel
from source.utils import plot_final_sizes, plot_and_save_all, plot_final_sizes_in_one


def main():
    # Solve controlled SEIR
    solve_controlled_seir()


def solve_controlled_seir():
    # Hubei & China population
    hubei_population = 58160000.0
    chn_population = 1437000000.0
    population = chn_population - hubei_population
    # Set time intervals
    incubation_period = 5.1
    # List of R0s
    r_0_list = [2.1, 2.6, 3.1]
    # r_0_list = [3.1]
    # List of infectious periods associated to R0s
    infectious_period_list = [1.7, 3.3, 5.6]
    # infectious_period_list = [5.6]
    # Initial values for exposed classes
    exposed_class_init = [[24163.0, 4619.0], [27504.0, 4538.0], [32677.0, 4088.0]]
    # exposed_class_init = [[32677.0, 4088.0]]
    # Initial values for infected classes
    infected_class_init = [[4619.0, 18.0, 14.0], [4538.0, 2.0, 1.0], [4088.0, 90.0, 0.0]]
    # infected_class_init = [[4088.0, 90.0, 0.0]]
    # Initial list for final sizes determined per R0 and t*
    final_sizes = []
    # List of t*-s
    t_stars = np.arange(20, 61, 1)
    for (r_0, infectious_period, exposed_init, infected_init) in \
            zip(r_0_list, infectious_period_list, exposed_class_init, infected_class_init):
        # Set rates
        beta = r_0 / (population * infectious_period)
        alpha = 1 / incubation_period
        gamma = 1 / infectious_period
        params = np.array((beta, alpha, gamma))
        # Set initial values
        init_values = {"e0": exposed_init,
                       "i0": infected_init,
                       "r0": [0]
                       }
        not_susceptible = sum([item for sublist in init_values.values() for item in sublist])
        init_values.update({"s0": [population - not_susceptible]})
        # Initial list for final sizes for all t_star
        fs = []
        for t_star in t_stars:
            # Instantiate epidemic model
            model = EpidemicModel(init_values, t_star=t_star, r_0=r_0)
            x0 = np.array(model.get_initial_values())
            # Solve epidemic model
            t = np.linspace(0, 200, 100000)
            solution = solve_model(t, x0, params, model)
            # solution = solve_model(t, x0, params, model)
            # Save final size for current t_star (recovered cases)
            # fs.append(solution.y[-1, -1])
            fs.append(solution[-1, -1])

            # Plot solutions and save figures
            # plot_and_save_all(t, solution.y.T, r_0, t_star)
            plot_and_save_all(t, solution, r_0, t_star)

        # Store final sizes for current R0
        final_sizes.append(fs)

    # Plot and save t_star vs. final size figures for all R0s
    for (r_0, final_size) in zip(r_0_list, final_sizes):
        plot_final_sizes(r_0, t_stars, final_size)
    plot_final_sizes_in_one(t_stars, final_sizes)   

    r_0_dict = {r_0_list[idx]: idx for idx in range(len(r_0_list))}
    return t_stars, final_sizes, r_0_dict


def solve_model(t, x0, params, model):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(model.get_model, x0, t, args=(params,))
    return x


def solve_model_2(t, x0, params, model):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = solve_ivp(fun=model.get_model_2, t_span=[min(t), max(t)], y0=x0, t_eval=t, args=(params,))
    return x


if __name__ == "__main__":
    main()
