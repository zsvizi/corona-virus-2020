import numpy as np
import itertools
from scipy import optimize
from scipy.stats import binom
from source.main_SEIR import solve_controlled_seir
from source.risk_test import test_function


def get_heatmap(max_number_summands):
    # Compute z
    r_locs = np.arange(1.0, 4.0, 0.01)
    extinction_probability = compute_extinction_probability(r_locs)

    # Get final sizes via solving SEIR models
    r0 = 2.6
    final_size, t_stars = get_final_size(r0)
    connectivity_ratio = np.arange(0, 0.2, 0.01)

    # Compute risks for possible combinations
    risk_ravel, order = compute_risk(final_sizes=final_size,
                                     connectivity_ratio=connectivity_ratio,
                                     extinction_probability=extinction_probability,
                                     max_number=max_number_summands)

    # Get all combinations of variables
    all_combinations = get_combinations(r_locs=r_locs,
                                        connectivity_ratio=connectivity_ratio,
                                        t_stars=t_stars, order=order)

    heatmap = {"heatmap": np.append(all_combinations, risk_ravel, axis=1), "t_stars": t_stars,
               "theta": connectivity_ratio, "r_locs": r_locs}
    return heatmap, order


def compute_risk(final_sizes, connectivity_ratio, extinction_probability, max_number=100):
    # DIMENSIONS: t_star x connectivity (theta) x r_loc (~z)
    risk = np.zeros((len(final_sizes), len(connectivity_ratio), len(extinction_probability)))
    j, k = [0, 0]
    for fs in final_sizes:
        k = 0
        for theta in connectivity_ratio:
            # Calculate binomial probabilities
            probs = binom.pmf(np.arange(max_number), fs, theta)
            # Calculate vector of risks corresponding to a specified fs and theta
            risk[j, k, :] = np.dot(np.array(probs), 1 - extinction_probability ** np.arange(max_number).reshape(-1, 1))
            k += 1
        j += 1
    # Order of previous loops determines order of variables - change IFF loop is changed!
    order = {"r_locs": 2, "connectivity_ratio": 1, "t_stars": 0}
    risk_ravel = risk.reshape((-1, 1))
    return risk_ravel, order


def get_final_size(r0):
    # Get t_stars, final_sizes, possible R0s
    t_stars, final_sizes, r_0_dict = solve_controlled_seir()
    # Get final sizes corresponding to the specified r0
    final_size_idx = r_0_dict[r0]
    final_size = final_sizes[final_size_idx]
    return final_size, t_stars


def compute_extinction_probability(r_locs):
    extinction_probability = \
        optimize.fixed_point(generator_function_neg_binom, 0.5 * np.ones(r_locs.shape[0]), args=(r_locs,))
    return extinction_probability


def generator_function_neg_binom(z, r_loc):
    # Parameters of negative binomial distribution: k (dispersion), p
    # Mean: (1-p)*k/p = r_loc -> p = k / (k + r_loc)
    # PGF:((1-p) / (1 - p*z))**k
    k = 2.2
    p = k / (k + r_loc)
    return np.power(p / (1 - (1 - p) * z), k)


def get_combinations(r_locs, t_stars, connectivity_ratio, order):
    # DIMENSIONS: t_star x connectivity (theta) x r_loc (~z)
    all_combinations_array = [[], [], []]
    all_combinations_array[order["t_stars"]] = t_stars
    all_combinations_array[order["connectivity_ratio"]] = connectivity_ratio
    all_combinations_array[order["r_locs"]] = r_locs
    all_combinations = np.array(list(itertools.product(*all_combinations_array)))
    return all_combinations


def main():
    max_number_summands = 100
    # Get heat map
    heatmap, order = get_heatmap(max_number_summands)
    heat_map = heatmap["heatmap"]
    t_stars = heatmap["t_stars"]
    connectivity_ratio = heatmap["theta"]
    r_locs = heatmap["r_locs"]
    # Check shape of heatmap
    print(heat_map.shape)

    test_function(heat_map=heat_map,
                  t_stars=t_stars,
                  connectivity_ratio=connectivity_ratio,
                  r_locs=r_locs,
                  max_number_summands=max_number_summands,
                  order=order)


if __name__ == "__main__":
    main()

