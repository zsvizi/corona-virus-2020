from __future__ import division
from lmfit import report_fit
import matplotlib.pyplot as plt
import numpy as np

from source.data_handling import load_data, generate_data
from source.fit import fit_model, solve_model
from source.model import R_0, number_of_E_compartments, number_of_I_compartments, EpidemicModel, \
    chn_population

plt.style.use('ggplot')


def main():
    bool_fit = False
    bool_solve_orig_seir = True
    bool_solve_controlled_seir = True

    if bool_fit:
        fitting()

    init_for_control = []
    if bool_solve_orig_seir:
        init_for_control = solve_orig_seir()

    if bool_solve_controlled_seir:
        solve_controlled_seir(init_for_control)


def solve_orig_seir():
    # hubei_population = 59170000
    # Set rates for generated data
    incubation_period = 5
    infectious_period = 10
    beta = R_0 / (chn_population * infectious_period)
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
    plot_cumulative(cumulative, t)
    plot_all(solution, t)

    init_for_control = print_endpoint(solution, cumulative, t)
    return init_for_control


def solve_controlled_seir(init_for_control):
    hubei_population = 59170000
    # Hopkins: 23.01.2020. total number of infected = 865, infected cases NOT from Hubei province = 316
    outside_hubei_ratio = 316 / 865
    init = outside_hubei_ratio * init_for_control
    s0 = chn_population - (hubei_population - np.sum((1 - outside_hubei_ratio) * init_for_control[1:-1]))
    # Set rates for generated data
    incubation_period = 5
    infectious_period = 10
    r_0_list = [1.95, 2.15, 2.13, 2.2, 2.6, 2.9]
    solutions = []
    for r_0 in r_0_list:
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
        # Instantiate epidemic model
        sols = []
        for t_star in np.arange(10, 51, 1):
            model = EpidemicModel(init_values, t_star=t_star, r_0=r_0)
            x0 = np.array(model.get_initial_values())
            r0_to_reach = 0.1
            endtime = t_star * (r_0 - r0_to_reach) / (r_0 - 1)
            # Solve epidemic model
            t = np.linspace(0, 10 * t_star, 100000)
            solution = solve_model(t, x0, params, model)
            sols.append(solution)
            plot_and_save_fig(t, solution, r_0, t_star)
        solutions.append(sols)
    cases_after_23jan = {"0": 316, "0.5": 367, "1": 591,
                         "1.5": 638, "2": 1004, "2.5": 1314, "3": 1402, "3.5": 1815}
    for idx in np.arange(len(r_0_list)):
        best_sol_idx = find_closest_solution(solutions[idx], cases_after_23jan)
        print("Best t* for R0=", r_0_list[idx], "is:", np.arange(10, 51, 1)[best_sol_idx])


def find_closest_solution(solution, cases: dict):
    best_sol_idx = 0
    t_star = 10
    for sol in solution:
        for t in cases.keys():
            index = (10 * t_star)/100000
        t_star += 1
    return best_sol_idx


def print_endpoint(solution, cumulative, t):
    condition = abs(cumulative - 830) < 10 ** (-1)
    element = cumulative[condition]
    element = np.min(abs(element - 830)) + 830
    index = list(cumulative).index(element)
    print("Time (days):", t[index])
    print("E1:", solution[index, 1])
    print("E2:", solution[index, 2])
    print("I1:", solution[index, 3])
    print("I2:", solution[index, 4])
    print("I3:", solution[index, 5])
    print("C:", solution[index, -1])
    return solution[index]


def plot_cumulative(cumulative, t):
    plt.plot(t, cumulative, '-', linewidth=2)
    plt.xlabel("Time (days)")
    plt.ylabel("Cumulative number of infected cases")
    plt.savefig("..\\data\\original_model\\solution_SEIR_orig.png", format="png")
    plt.savefig("..\\data\\original_model\\solution_SEIR_orig.eps", format="eps")
    # plt.show()


def plot_all(solution, t):
    plt.plot(t, solution[:, 1:], linewidth=2)
    plt.xlabel("Time (days)")
    plt.ylabel("State values")
    plt.legend(('E1(t)', 'E2(t)', 'I1(t)', 'I2(t)', 'I3(t)', 'R(t)', 'C(t)'), loc='upper left')
    plt.savefig("..\\data\\original_model\\solution_SEIR_orig_all.png", format="png")
    plt.savefig("..\\data\\original_model\\solution_SEIR_orig_all.eps", format="eps")
    # plt.show()


def plot_and_save_fig(t, solution, r_0, t_star):
    plt.plot(t, solution[:, 1:], linewidth=2)
    plt.xlabel("Time (days)")
    plt.ylabel("State values")
    plt.legend(('E1(t)', 'E2(t)', 'I1(t)', 'I2(t)', 'I3(t)', 'R(t)', 'C(t)'), loc='upper left')
    plt.savefig("..\\data\\controlled_models\\solution_SEIR_control_R0_" + str(r_0) + "_t_" + str(t_star) + "_all.png",
                format="png")
    plt.savefig("..\\data\\controlled_models\\solution_SEIR_control_R0_" + str(r_0) + "_t_" + str(t_star) + "_all.eps",
                format="eps")
    plt.close()
    print("R0 =", r_0, "t_star =", t_star)


def fitting():
    # Get data
    path = '..\\data\\cases_first_period.txt'
    t, data = get_data(path)
    # Get data to fit
    data_to_fit = np.diff(data[:, 1:2], axis=0)
    init_values = get_init_values(data, data_to_fit)
    # fit to generated data
    final, result = fit_model(data, init_values, t)
    # plot data and fitted curves
    plt.plot(t, data, 'o')
    plt.plot(t, final, '-', linewidth=2)
    plt.savefig("..\\data\\fitting.png", format="png")
    plt.show()
    # display fitted statistics
    report_fit(result)


def get_data(path):
    # Load data
    t, data = load_data(path)
    test = False
    if test:
        # generate data
        t, data = generate_data()
    return t, data


def get_init_values(data, data_to_fit):
    init_values = {"e0": [0, 0],
                   "i0": [0, data[0, 1]],
                   "r0": [data[0, 2]]}
    init_value_s = {"s0": [chn_population -
                           (sum(init_values["e0"]) + sum(init_values["i0"]) + sum(init_values["r0"]))],
                    "c0": [data_to_fit[0, 0]]}
    init_values.update(init_value_s)
    return init_values


if __name__ == "__main__":
    main()
