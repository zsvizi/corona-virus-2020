from scipy import optimize
from scipy.stats import binom
import math
import random
import numpy as np
import itertools
import sys

'''
risk.py calculates the risk of an epidemic in a destination country.

'''

def get_heatmap(c, theta, r_loc):
    r_stars = c
    connectivity_ratio = theta
    r_locs = r_loc
    
    z = compute_z(r_locs)
    
    # Compute risks for possible combinations
    risk_ravel, order = compute_risk(r_stars=r_stars, connectivity_ratio=connectivity_ratio, solution=z)

    # Get all combinations of variables
    all_combinations = get_combinations(r_locs=r_locs,
                                        connectivity_ratio=connectivity_ratio,
                                        r_stars=r_stars, order=order)

    heatmap = {"heatmap": np.append(all_combinations, risk_ravel, axis=1), "r_stars": r_stars,
               "theta": connectivity_ratio, "r_locs": r_locs}
    return heatmap, order


def compute_z(r_locs):
    solution = optimize.fixed_point(generator_function_neg_binom, 0.5 * np.ones(r_locs.shape[0]), args=(r_locs,))
    return solution


def generator_function_neg_binom(z, r_loc):
    # Parameters: k (dispersion), p
    # Mean: (1-p)*k/p = r_loc -> p = k / (k + r_loc)
    # PGF:((1-p) / (1 - p*z))**k
    # Dispersion parameter - TO BE DETERMINED
    k = 0.64
    p = k / (k + r_loc)
    return np.power(p / (1 - (1 - p) * z), k)

def compute_risk(r_stars, connectivity_ratio, solution):
    # DIMENSIONS: r_star (final size at t_star, later maybe t_star) x connectivity (theta) x r_loc (~z)
    risk = np.zeros((len(r_stars), len(connectivity_ratio), len(solution)))
    allin = len(r_stars)*len(connectivity_ratio)
    k=0
    for (r_starIndex, r_star) in enumerate(r_stars):
        for (thetaIndex, theta) in enumerate(connectivity_ratio):
            sys.stdout.write("\033[F")
            sys.stdout.write("\033[K")
            print(str(allin) + "/" + str(k))
            risk[r_starIndex, thetaIndex, :] = 1 - ((1 - theta) + theta * solution) ** r_star
            k+=1
    # Order of previous loops determines order of variables
    order = {"r_locs": 2, "connectivity_ratio": 1, "r_stars": 0}
    risk_ravel = risk.reshape((-1, 1))
    return risk_ravel, order


def get_combinations(r_locs, r_stars, connectivity_ratio, order):
    # DIMENSIONS: r_star (final size at t_star, later maybe t_star) x connectivity (theta) x r_loc (~z)
    all_combinations_array = [[], [], []]
    all_combinations_array[order["r_stars"]] = r_stars
    all_combinations_array[order["connectivity_ratio"]] = connectivity_ratio
    all_combinations_array[order["r_locs"]] = r_locs
    all_combinations = np.array(list(itertools.product(*all_combinations_array)))
    return all_combinations


def write_file(connectivity_ratio, heat_map, r_locs, r_stars):
    with open("../data/heatmap_to_save.txt", "a+") as file:
        file.write(str(r_locs)[1:-1])
        file.write("\n")

        file.write(str(connectivity_ratio)[1:-1])
        file.write("\n")

        file.write(str(r_stars)[1:-1])
        file.write("\n")

        for vector in heat_map:
            file.write(str(vector)[1:-1])
            file.write("\n")

            
def main():
    # Get heat map
    heatmap, order = get_heatmap()
    heat_map = heatmap["heatmap"]
    r_stars = heatmap["r_stars"]
    connectivity_ratio = heatmap["theta"]
    r_locs = heatmap["r_locs"]

if __name__ == "__main__":
    main()
