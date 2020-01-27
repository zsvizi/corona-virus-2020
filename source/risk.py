from scipy import optimize
import numpy as np
import itertools
"""
Risk = 1 - z ** p

For getting z:
    Galton-Watson process:
    - generator function of offspring distribution: g(z)
    - z is the solution of equation z = g(z)
    => for Poisson distribution: z = e ** (R_loc * (z - 1))

For getting p (= number of branches):
    - p = (1-s) * R(t_*) * theta
    - p depends on
        - t_*
        - theta
        
=> Risk = Risk(R_loc, t_*, theta)

Goal: plot heat maps, where heat value = Risk for 
    - R_loc vs. t_* and fixed theta
    - R_loc vs. theta and fixed t_*
    - t_* vs. theta and fixed R_loc
    
Extremal cases:
    > R_loc -> 1
    > theta -> 0
"""


def generator_function(z, r_loc):
    return np.exp(r_loc * (z - 1))


def compute_risk(number_of_branches, solution):
    risk_arrays = [1 - x ** number_of_branches for x in solution]
    risk = np.stack(risk_arrays, axis=0)
    risk_ravel = risk.reshape((-1, 1))
    return risk_ravel


def get_combinations(connectivity_ratio, r_locs, r_stars):
    all_combinations_array = [r_locs, r_stars, connectivity_ratio]
    all_combinations = np.array(list(itertools.product(*all_combinations_array)))
    return all_combinations


def compute_p(r_stars, connectivity_ratio):
    screening_efficacy = 0.9
    number_of_branches = (1 - screening_efficacy) * np.outer(r_stars, connectivity_ratio)
    return number_of_branches


def compute_z(r_locs):
    solution = optimize.fixed_point(generator_function, 0.5 * np.ones(r_locs.shape[0]), args=(r_locs,))
    return solution


def main():
    # Compute z
    r_locs = np.array([0.8, 1.2, 1.6, 2.0])
    z = compute_z(r_locs)

    # Compute p
    r_stars = np.array([1000, 1500, 2500])
    connectivity_ratio = np.array([0.2, 0.5, 0.75, 0.9, 0.95])
    p = compute_p(r_stars, connectivity_ratio)

    # Get all combinations of variables
    all_combinations = get_combinations(connectivity_ratio, r_locs, r_stars)

    # Compute risks for possible combinations
    risk_ravel = compute_risk(p, z)

    # Get heat map
    heat_map = np.append(all_combinations, risk_ravel, axis=1)
    print(heat_map)


if __name__ == "__main__":
    main()

