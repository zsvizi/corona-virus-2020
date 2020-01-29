from lmfit import Parameters
import copy
number_of_S_compartments = 1
number_of_E_compartments = 2
number_of_I_compartments = 3
number_of_R_compartments = 1
chn_population = 1436980159
R_0 = 2.6


class EpidemicModel:
    def __init__(self, init_values=None, t_star=None, r_0=None):
        if init_values is None:
            init_values = {"s0": [100], "e0": [1, 0], "i0": [0, 0, 0], "r0": [0], "c0": [0]}
        self.susceptible = Susceptible(init_values["s0"])
        self.exposed = Exposed(init_values["e0"])
        self.infected = Infected(init_values["i0"])
        self.recovered = Recovered(init_values["r0"])
        self.cumulative = Cumulative(init_values["c0"])
        self.params = Parameters()
        self.parameter_init()
        self.t_star = t_star
        self.r_0 = r_0

    def get_model(self, xs, t, ps):
        """
        Epidemic model
        """
        try:
            alpha, gamma = ps['alpha'].value, ps['gamma'].value
            beta = ps['beta'].value
        except:
            beta, alpha, gamma = ps
        # beta = R_0 * gamma / chn_population

        s, e1, e2, i1, i2, i3, r, c = xs
        if self.t_star is None:
            control = 1
        else:
            r0_to_reach = 0.5
            slope = (self.r_0 - 1) / (self.r_0 * self.t_star)
            y_to_reach = (self.r_0 - r0_to_reach) / self.r_0
            control = 1 - min(y_to_reach, slope * t)
        model_eq = [
            -beta * control * s * (i1 + i2 + i3),  # S'(t)
            beta * control * s * (i1 + i2 + i3) - alpha * e1,  # E1'(t)
            alpha * e1 - alpha * e2,  # E2'(t)
            alpha * e2 - gamma * i1,  # I1'(t)
            gamma * i1 - gamma * i2,  # I2'(t)
            gamma * i2 - gamma * i3,  # I3'(t)
            gamma * i3,  # R(t)
            beta * s * (i1 + i2 + i3)  # C'(t)
        ]
        return model_eq

    def add_rates(self):
        self.params.add('alpha', value=1 / (5 * number_of_E_compartments),
                        min=1 / (10 * number_of_E_compartments),
                        max=1 / (3 * number_of_E_compartments))
        self.params.add('gamma', value=1 / (6 * number_of_I_compartments),
                        min=1 / (9 * number_of_I_compartments),
                        max=1 / (3 * number_of_I_compartments))
        self.params.add('beta', value=2.6 * 5 * self.params['gamma']/chn_population,
                        min=1.4 * 3 * self.params['gamma']/chn_population,
                        max=3 * 9 * self.params['gamma']/chn_population)

    def get_initial_values(self):
        init_values = copy.deepcopy(self.susceptible.get_initial_values())
        init_values.extend(self.exposed.get_initial_values())
        init_values.extend(self.infected.get_initial_values())
        init_values.extend(self.recovered.get_initial_values())
        init_values.extend(self.cumulative.get_initial_values())
        return init_values

    def parameter_init(self):
        self.add_initial_values()
        self.add_rates()

    def add_initial_values(self):
        self.susceptible.add_init_value(self.params)
        self.exposed.add_init_value(self.params)
        self.infected.add_init_value(self.params)
        self.recovered.add_init_value(self.params)
        self.cumulative.add_init_value(self.params)


class Susceptible:
    def __init__(self, s0):
        self.s0 = s0

    def get_initial_values(self):
        return self.s0

    def add_init_value(self, ps: Parameters):
        for i in range(len(self.s0)):
            ps.add('s' + str(i+1) + '0', value=float(self.s0[i]), vary=False)


class Exposed:
    def __init__(self, e0):
        self.e0 = e0

    def get_initial_values(self):
        return self.e0

    def add_init_value(self, ps: Parameters):
        for i in range(number_of_E_compartments):
            ps.add('e' + str(i+1) + '0', value=float(self.e0[i]), vary=False)


class Infected:
    def __init__(self, i0):
        self.i0 = i0

    def get_initial_values(self):
        return self.i0

    def add_init_value(self, ps: Parameters):
        for i in range(number_of_I_compartments):
            # ps.add('i' + str(i+1) + '0', value=float(self.i0[i]), min=0, max=200)
            ps.add('i' + str(i + 1) + '0', value=float(self.i0[i]), vary=False)


class Recovered:
    def __init__(self, r0):
        self.r0 = r0

    def get_initial_values(self):
        return self.r0

    def add_init_value(self, ps: Parameters):
        for i in range(len(self.r0)):
            ps.add('r' + str(i+1) + '0', value=float(self.r0[i]), vary=False)


class Cumulative:
    def __init__(self, c0):
        self.c0 = c0

    def get_initial_values(self):
        return self.c0

    def add_init_value(self, ps: Parameters):
        for i in range(len(self.c0)):
            ps.add('c' + str(i+1) + '0', value=float(self.c0[i]), vary=False)
