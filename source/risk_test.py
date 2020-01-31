import numpy as np
import random
from scipy.stats import binom


def test_function(heat_map, t_stars, connectivity_ratio, r_locs, max_number_summands, order):
    passed = True
    number_of_tests = 1000
    r_0 = 2.6

    # Due to circular dependency
    from source.risk import get_final_size, compute_extinction_probability
    final_size, _ = get_final_size(r_0)

    for _ in np.arange(number_of_tests):
        cr_len = len(connectivity_ratio)
        rl_len = len(r_locs)
        fs_len = len(t_stars)
        lengths = [0, 0, 0]
        lengths[order["t_stars"]] = fs_len
        lengths[order["connectivity_ratio"]] = cr_len
        lengths[order["r_locs"]] = rl_len

        cr_index = random.randint(0, cr_len-1)
        rl_index = random.randint(0, rl_len-1)
        fs_index = random.randint(0, fs_len-1)
        indices = [0, 0, 0]
        indices[order["t_stars"]] = fs_index
        indices[order["connectivity_ratio"]] = cr_index
        indices[order["r_locs"]] = rl_index

        cr_element = connectivity_ratio[cr_index]
        rl_element = r_locs[rl_index]
        fs_element = final_size[fs_index]

        z = compute_extinction_probability(np.array([rl_element]))
        expected = np.dot(np.array(binom.pmf(np.arange(max_number_summands), fs_element, cr_element)),
                          1 - z ** np.arange(max_number_summands))
        res_index = ((fs_index * cr_len * rl_len) + (cr_index * rl_len) + rl_index)
        result = heat_map[res_index, -1]
        if abs(result - expected) > 10**(-5):
            passed = False
            print("Test failed")
            break
    if passed:
        print("Test passed!")