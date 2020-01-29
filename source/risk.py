from scipy import optimize
from scipy.stats import binom
import math
import numpy as np
import itertools


def get_heatmap(max_number_summands):
    # Compute z
    r_locs = np.arange(1.0, 4.0, 0.01)
    z = compute_z(r_locs)

    # Compute p
    # No r(t_star) data yet
    r_stars = np.arange(1000, 10000, 100)
    connectivity_ratio = np.arange(0, 0.2, 0.01)

    # Get all combinations of variables
    all_combinations = get_combinations(r_locs=r_locs,
                                        connectivity_ratio=connectivity_ratio,
                                        r_stars=r_stars)

    # Compute risks for possible combinations
    risk_ravel = compute_risk(r_stars=r_stars, connectivity_ratio=connectivity_ratio, solution=z,
                              max_number=max_number_summands)

    return np.append(all_combinations, risk_ravel, axis=1), connectivity_ratio, r_locs, r_stars


def compute_z(r_locs):
    solution = optimize.fixed_point(generator_function_neg_binom, 0.5 * np.ones(r_locs.shape[0]), args=(r_locs,))
    return solution


def generator_function_neg_binom(z, r_loc):
    # Parameters: k (dispersion), p
    # Mean: p*k/(1-p) = r_loc -> p = r_loc / (k + r_loc)
    # PGF:((1-p) / (1 - p*z))**k
    # Dispersion parameter - TO BE DETERMINED
    k = 2.2
    p = r_loc / (k + r_loc)
    return np.power((1 - p) / (1 - p * z), k)


def compute_risk(r_stars, connectivity_ratio, solution, max_number=100):
    # DIMENSIONS: r_stars x theta x z
    risk = np.zeros((len(r_stars), len(connectivity_ratio), len(solution)))
    j, k = [0, 0]
    for r_star in r_stars:
        k = 0
        for theta in connectivity_ratio:
            probs = binom.pmf(np.arange(max_number), r_star, theta)
            risk[j, k, :] = np.dot(np.array(probs), 1 - solution ** np.arange(max_number).reshape(-1, 1))
            k += 1
        j += 1
    # risk_arrays = [1 - x ** r_stars for x in z]
    # risk = np.stack(risk_arrays, axis=0)
    risk_ravel = risk.reshape((-1, 1))
    return risk_ravel


def get_combinations(r_locs, r_stars, connectivity_ratio):
    all_combinations_array = [r_stars, connectivity_ratio, r_locs]
    all_combinations = np.array(list(itertools.product(*all_combinations_array)))
    return all_combinations


def write_file(connectivity_ratio, heat_map, r_locs, r_stars):
    with open("..\\data\\teszt_2.txt", "a+") as file:
        file.write(str(r_locs)[1:-1])
        file.write("\n")

        file.write(str(connectivity_ratio)[1:-1])
        file.write("\n")

        file.write(str(r_stars)[1:-1])
        file.write("\n")

        for vector in heat_map:
            file.write(str(vector)[1:-1])
            file.write("\n")


def generator_function_poisson(z, r_loc):
    return np.exp(r_loc * (z - 1))


def compute_p(r_stars, connectivity_ratio):
    screening_efficacy = 0.9
    number_of_branches = (1 - screening_efficacy) * np.outer(r_stars, connectivity_ratio)
    return number_of_branches


def test(heat_map, connectivity_ratio, r_locs, r_stars, max_number_summands):
    cr_index = int(math.ceil(len(connectivity_ratio)/2))
    rl_index = int(math.ceil(len(r_locs) / 2))
    rs_index = int(math.ceil(len(r_stars) / 2))
    cr_element = connectivity_ratio[cr_index]
    rl_element = r_locs[rl_index]
    rs_element = r_stars[rs_index]
    z = compute_z(np.array([rl_element]))
    expected = np.dot(np.array(binom.pmf(np.arange(max_number_summands), rs_element, cr_element)),
                      1 - z ** np.arange(max_number_summands))
    res_index = ((rs_index * len(connectivity_ratio) * len(r_locs)) + (cr_index * len(r_locs)) + rl_index)
    result = heat_map[res_index, -1]
    print("Expected:", expected)
    print("Result:", result)


def main():
    max_number_summands = 100
    # Get heat map
    heat_map, connectivity_ratio, r_locs, r_stars = get_heatmap(max_number_summands)
    print(heat_map.shape)

    test(heat_map, connectivity_ratio, r_locs, r_stars, max_number_summands)
    # write_file(connectivity_ratio, heat_map, r_locs, r_stars)


if __name__ == "__main__":
    main()

