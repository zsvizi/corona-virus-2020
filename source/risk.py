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


def generator_function_poisson(z, r_loc):
    return np.exp(r_loc * (z - 1))


def generator_function_neg_binom(z, r_loc):
    # (p*z / (1 - (1-p)*z))**n
    # p*n / (1-p) = r_loc
    p = 0.5
    return np.power((p*z) / (1 - (1-p)*z), (r_loc * (1-p)) / p)


def compute_risk(number_of_branches, solution):
    risk_arrays = [1 - x ** number_of_branches for x in solution]
    risk = np.stack(risk_arrays, axis=0)
    risk_ravel = risk.reshape((-1, 1))
    return risk_ravel


def get_combinations(r_locs, r_stars, connectivity_ratio):
    all_combinations_array = [r_locs, r_stars, connectivity_ratio]
    all_combinations = np.array(list(itertools.product(*all_combinations_array)))
    return all_combinations


def compute_p(r_stars, connectivity_ratio):
    screening_efficacy = 0.9
    number_of_branches = (1 - screening_efficacy) * np.outer(r_stars, connectivity_ratio)
    return number_of_branches


def compute_z(r_locs):
    solution = optimize.fixed_point(generator_function_poisson, 0.5 * np.ones(r_locs.shape[0]), args=(r_locs,))
    return solution


def get_heatmap():
    # Compute z
    r_locs = np.arange(1.05, 1.09, 0.05)
    z = compute_z(r_locs)

    # Compute p
    # t_stars = np.arange(10, 40, 2)
    r_stars = np.arange(1000, 2000, 1)
    connectivity_ratio = np.arange(0, 0.2, 0.01)
    p = compute_p(r_stars, connectivity_ratio)

    # Get all combinations of variables
    # Dimensions: r_locs x connectivity_ratio x r_stars
    all_combinations = get_combinations(r_locs=r_locs,
                                        connectivity_ratio=connectivity_ratio,
                                        r_stars=r_stars)

    # Compute risks for possible combinations
    risk_ravel = compute_risk(p, z)

    return np.append(all_combinations, risk_ravel, axis=1), r_locs, connectivity_ratio, r_stars


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


def main():

    # Get heat map
    heat_map, connectivity_ratio, r_locs, r_stars = get_heatmap()
    print(heat_map.shape)

    write_file(connectivity_ratio, heat_map, r_locs, r_stars)


if __name__ == "__main__":
    main()

