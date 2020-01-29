import numpy as np
from lmfit import minimize
from scipy.integrate import odeint

from source.model import EpidemicModel


def fit_model(data, init_values, t):
    # Instantiate epidemic model
    model = EpidemicModel(init_values)
    # fit model and find predicted values
    final, result = fit_and_predict(data, t, model)
    return final, result


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