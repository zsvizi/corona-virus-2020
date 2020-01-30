import csv
import numpy as np

from source.main_SEIR import solve_model
from source.model import number_of_E_compartments, number_of_I_compartments, EpidemicModel


def get_data(path):
    # Load data
    t, data = load_data(path)
    test = False
    if test:
        # generate data
        t, data = generate_data()
    return t, data


def generate_data():
    # Set rates for generated data
    beta, alpha, gamma = 0.7, 0.25/number_of_E_compartments, 0.5/number_of_I_compartments
    true_params = np.array((beta, alpha, gamma))
    # Instantiate epidemic model
    model = EpidemicModel()
    # Get initial values
    x0 = np.array(model.get_initial_values())
    # Solve epidemic model
    t = np.linspace(0, 10, 10)
    data = solve_model(t, x0, true_params, model)
    # Perturbate with random noise
    data += np.random.normal(size=data.shape)
    return t, data[:, -2:]


def load_data(path):
    with open(path, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        data = np.array(list(reader)).astype(float)
    t = data[1:, 0]
    # data_to_fit = data[:, 1:2]
    return t, data
