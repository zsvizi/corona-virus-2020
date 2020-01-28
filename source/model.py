from lmfit import Parameters


class Susceptible:
    def __init__(self, s0):
        self.s0 = s0

    def get_initial_values(self):
        return self.s0

    def add_init_value(self, ps: Parameters):
        ps.add('s0', value=float(self.s0), vary=False)


class Exposed:
    def __init__(self, e0):
        self.e0 = e0

    def get_initial_values(self):
        return self.e0

    def add_init_value(self, ps: Parameters):
        ps.add('e0', value=float(self.e0), vary=False)


class Infected:
    def __init__(self, i0):
        self.i0 = i0

    def get_initial_values(self):
        return self.i0

    def add_init_value(self, ps: Parameters):
        ps.add('i0', value=float(self.i0), vary=False)


class Recovered:
    def __init__(self, r0):
        self.r0 = r0

    def get_initial_values(self):
        return self.r0

    def add_init_value(self, ps: Parameters):
        ps.add('r0', value=float(self.r0), vary=False)


class EpidemicModel:
    def __init__(self, init_values=None):
        if init_values is None:
            init_values = {"s0": 100, "e0": 1, "i0": 0, "r0": 0}
        self.susceptible = Susceptible(init_values["s0"])
        self.exposed = Exposed(init_values["e0"])
        self.infected = Infected(init_values["i0"])
        self.recovered = Recovered(init_values["r0"])
        self.params = Parameters()
        self.parameter_init()

    @staticmethod
    def get_model(xs, t, ps):
        """
        Epidemic model
        """
        t_star = 5
        try:
            beta, alpha, gamma = ps['beta'].value, ps['alpha'].value, ps['gamma'].value
        except:
            beta, alpha, gamma = ps
        s, e, i, r = xs
        return [
            -beta * (1 - min(0.61, 0.61 / t_star * t)) * s * i,
            beta * (1 - min(0.61, 0.61 / t_star * t)) * s * i - alpha * e,
            alpha * e - gamma * i,
            gamma * i
        ]

    def get_initial_values(self):
        return self.susceptible.get_initial_values(), \
               self.exposed.get_initial_values(), \
               self.infected.get_initial_values(), \
               self.recovered.get_initial_values(),

    def parameter_init(self):
        self.add_initial_values()
        self.add_rates()

    def add_initial_values(self):
        self.susceptible.add_init_value(self.params)
        self.exposed.add_init_value(self.params)
        self.infected.add_init_value(self.params)
        self.recovered.add_init_value(self.params)

    def add_rates(self):
        self.params.add('beta', value=0.2, min=0, max=3)
        self.params.add('alpha', value=0.05, min=0, max=1)
        self.params.add('gamma', value=0.05, min=0, max=1)