import math

import numpy as np
from matplotlib import pyplot as plt


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


def plot_final_sizes(r_0, t_star, final_size):
    plt.plot(t_star, final_size, linewidth=2)
    plt.xlabel("t*")
    plt.ylabel("Final size")
    plt.savefig("..\\data\\final_size\\final_size_R0_" + str(r_0) + ".png",
                format="png")
    plt.savefig("..\\data\\final_size\\final_size_R0_" + str(r_0) + ".eps",
                format="eps")
    plt.close()
