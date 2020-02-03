from matplotlib import pyplot as plt
import numpy as np
plt.style.use('seaborn-whitegrid')


def plot_and_save_all(t, solution, r_0, t_star):
    plt.plot(t, solution[:, 1:], linewidth=2)
    plt.xlabel("Time (days)")
    plt.ylabel("State values")
    plt.legend(('E1(t)', 'E2(t)', 'I1(t)', 'I2(t)', 'I3(t)', 'R(t)'), loc='upper right')
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


def plot_final_sizes_in_one(t_stars, final_sizes):
    linestyles = [':', '-', '--']
    colors = ['g', 'blue', 'tomato']
    for (final_size, c, ls) in zip(final_sizes, colors, linestyles):
        plt.plot(t_stars, final_size, linewidth=2, color=c, linestyle=ls)
    plt.xlabel("$t_*$", fontsize=14)
    plt.ylabel("Final size", fontsize=14)
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MultipleLocator(2000000))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_function))
    plt.text(48, 10000000, r'$R_0 = 3.1$', fontsize=14, color='tomato')
    plt.text(53, 3800000, r'$R_0 = 2.6$', fontsize=14, color='blue')
    plt.text(53, 1500000, r'$R_0 = 2.1$', fontsize=14, color='g')
    plt.tight_layout()
    plt.savefig("..\\data\\final_size\\final_sizes_all.png",
                format="png", dpi=500)
    plt.savefig("..\\data\\final_size\\final_sizes_all.eps",
                format="eps", dpi=500)
    plt.close()


def format_function(value, tick_number):
    # find number of multiples of pi/2
    x = np.around(value / 10000000, decimals=2)
    return r"${0} \times 10^7$".format(x)

