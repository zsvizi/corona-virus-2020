import math

import numpy as np
from matplotlib import pyplot as plt


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
    plt.close()
    # plt.show()


def plot_and_save_all(t, solution, r_0, t_star):
    plt.plot(t, solution[:, 1:], linewidth=2)
    plt.xlabel("Time (days)")
    plt.ylabel("State values")
    plt.legend(('E1(t)', 'E2(t)', 'I1(t)', 'I2(t)', 'I3(t)', 'R(t)', 'C(t)'), loc='upper right')
    plt.savefig("..\\data\\controlled_models\\solution_SEIR_control_R0_" + str(r_0) + "_t_" + str(t_star) + "_all.png",
                format="png")
    plt.savefig("..\\data\\controlled_models\\solution_SEIR_control_R0_" + str(r_0) + "_t_" + str(t_star) + "_all.eps",
                format="eps")
    plt.close()
    print("R0 =", r_0, "t_star =", t_star)


def plot_and_save_compare(cases_data, r_0, solution, t, t_data, t_star):
    plt.plot(t_data, cases_data, 'o')
    t_max = np.max(t_data)
    t_max_idx = math.ceil(100000 * t_max / (10 * t_star))
    plt.plot(t[:t_max_idx + 1], solution[:t_max_idx + 1, 1:], linewidth=2)
    plt.xlabel("Time (days)")
    plt.ylabel("State values")
    plt.legend(('reported cases', 'E1(t)', 'E2(t)', 'I1(t)', 'I2(t)', 'I3(t)', 'R(t)', 'C(t)'), loc='upper right')
    plt.savefig("..\\data\\controlled_models_compare\\solution_SEIR_control_R0_" + str(r_0) + "_t_"
                + str(t_star) + "_compare.png", format="png")
    plt.savefig("..\\data\\controlled_models_compare\\solution_SEIR_control_R0_" + str(r_0) + "_t_"
                + str(t_star) + "_compare.eps", format="eps")
    plt.close()


def plot_final_sizes(r_0, final_size):
    plt.plot(np.arange(10, 51, 1), final_size, linewidth=2)
    plt.xlabel("t*")
    plt.ylabel("Final size")
    plt.savefig("..\\data\\final_size\\final_size_R0_" + str(r_0) + ".png",
                format="png")
    plt.savefig("..\\data\\final_size\\final_size_R0_" + str(r_0) + ".eps",
                format="eps")
    plt.close()