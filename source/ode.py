from __future__ import division
from lmfit import minimize, Parameters, report_fit
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')


def model_seir(xs, t, ps):
    """
    Epidemic model
    """
    try:
        beta = ps['beta'].value
        alpha = ps['alpha'].value
        gamma = ps['gamma'].value
    except:
        beta, alpha, gamma = ps

    s, e, i, r = xs
    return [-beta*s*i, beta*s*i - alpha*e, alpha*e - gamma*i, gamma*i]


def solution_seir(t, x0, ps):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(model_seir, x0, t, args=(ps,))
    return x


def residual_seir(ps, ts, data):
    x0 = ps['s0'].value, ps['e0'].value, ps['i0'].value, ps['r0'].value
    model = solution_seir(ts, x0, ps)
    return (model - data).ravel()


def seir_fit():
    t = np.linspace(0, 50, 100)
    x0 = np.array([100, 1, 0, 0])
    beta, alpha, gamma = 0.8, 0.1, 0.1
    true_params = np.array((beta, alpha, gamma))
    data = solution_seir(t, x0, true_params)
    data += np.random.normal(size=data.shape)
    # set parameters including bounds
    params = Parameters()
    params.add('s0', value=float(x0[0]), vary=False)
    params.add('e0', value=float(x0[1]), vary=False)
    params.add('i0', value=float(x0[2]), vary=False)
    params.add('r0', value=float(x0[3]), vary=False)
    params.add('beta', value=0.2, min=0, max=3)
    params.add('alpha', value=0.05, min=0, max=1)
    params.add('gamma', value=0.05, min=0, max=1)
    # fit model and find predicted values
    result = minimize(residual_seir, params, args=(t, data), method='leastsq')
    final = data + result.residual.reshape(data.shape)
    return data, final, result, t


def main():
    data, final, result, t = seir_fit()

    # plot data and fitted curves
    plt.plot(t, data, 'o')
    plt.plot(t, final, '-', linewidth=2)
    plt.show()

    # display fitted statistics
    report_fit(result)


if __name__ == "__main__":
    main()
