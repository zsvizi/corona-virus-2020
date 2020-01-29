from __future__ import division
from lmfit import minimize, report_fit
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import csv
from source.model import EpidemicModel, number_of_I_compartments, number_of_E_compartments, get_init_values

plt.style.use('ggplot')


def main():
    # fit to generated data
    data, final, result, t = fit_model()

    # plot data and fitted curves
    plot_result(data, final, t)

    # display fitted statistics
    report_fit(result)


def fit_model():
    # Load data
    t, data, init_values = load_data('..\\data\\cases_first_period.txt')
    test = False
    if test:
        # generate data
        t, data, init_values = generate_data()
    # Instantiate epidemic model
    model = EpidemicModel(init_values)
    # fit model and find predicted values
    final, result = fit_and_predict(data, t, model)
    return data, final, result, t


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
    return t, data[:, -2:], None


def load_data(path):
    with open(path, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        data = np.array(list(reader)).astype(float)
    t = data[1:, 0]
    init_values, data_to_fit = get_init_values(data)
    return t, data_to_fit, init_values


def fit_and_predict(data, t, model):
    # Execute fitting
    result = minimize(compute_residual, model.params, args=(t, data, model), method='leastsq')
    # Compute function values
    final = data + result.residual.reshape(data.shape)
    return final, result


def compute_residual(params, ts, data, model):
    x0 = model.get_initial_values()
    sol = solve_model(ts, x0, params, model)
    sol_fit = np.vstack([sol[0, -1:].reshape(-1, 1), np.diff(sol[:, -1], axis=0).reshape(-1, 1)])
    return (sol_fit - data).ravel()


def solve_model(t, x0, params, model):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(model.get_model, x0, t, args=(params,))
    return x


def plot_result(data, final, t):
    plt.plot(t, data, 'o')
    plt.plot(t, final, '-', linewidth=2)
    plt.savefig("..\\data\\fitting.png", format="png")
    plt.show()


if __name__ == "__main__":
    main()
