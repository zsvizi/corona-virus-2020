# To run this code, create a Bokeh server via
# `bokeh serve --show application`
# from folder ./corona-virus-2020

import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, RadioButtonGroup, Tabs, Panel, Select
from bokeh.plotting import figure

from .utils import EpidemicModel, solve_model


def model_solution(t, r_0, comp, population, incubation_period):
    # Set time intervals
    infectious_period = 1.7
    exposed_init = [24163.0, 4619.0]
    infected_init = [4619.0, 18.0, 14.0]
    t_star = 21

    beta = r_0 / (population * infectious_period)
    alpha = 1 / incubation_period
    gamma = 1 / infectious_period
    params = np.array((beta, alpha, gamma))
    # Set initial values
    init_values = {"e0": exposed_init,
                   "i0": infected_init,
                   "r0": [0]
                   }
    not_susceptible = sum([item for sublist in init_values.values() for item in sublist])
    init_values.update({"s0": [population - not_susceptible]})

    model = EpidemicModel(init_values, t_star=t_star, r_0=r_0)
    x0 = np.array(model.get_initial_values())
    # Solve epidemic model
    solution = solve_model(t, x0, params, model)
    sol = []
    if comp == 'R':
        sol = solution[:, -1]
    elif comp == 'I':
        sol = solution[:, -2]
    return sol


# Set up global variables
hubei_population = 58160000.0
chn_population = 1437000000.0
outside_china = chn_population - hubei_population
T = 200
pop_dict = {"Szeged": outside_china / 10000, "Budapest": outside_china, "Debrecen":outside_china / 500}
r0_list = [2.1, 2.8]

# Set up widgets
ip = Slider(title="Incubation period", value=5.1, start=5.0, end=14.0, step=0.1)
button = RadioButtonGroup(labels=["R0=2.1", "R0=2.8"], active=0)
select = Select(title="City:", value="Budapest", options=["Szeged", "Budapest", "Debrecen"])

# Initial plot
x = np.linspace(0, T, 100000)
y = model_solution(x,
                   comp='R',
                   r_0=r0_list[button.active],
                   incubation_period=ip.value,
                   population=pop_dict[select.value])
source_1 = ColumnDataSource(data=dict(x=x, y=y))
y = model_solution(x,
                   comp='I',
                   r_0=r0_list[button.active],
                   incubation_period=ip.value,
                   population=pop_dict[select.value])
source_2 = ColumnDataSource(data=dict(x=x, y=y))


# Set up plot
plot_1 = figure(plot_height=400, plot_width=400, title="SEIR model ",
                tools="crosshair,pan,reset,save,wheel_zoom",
                x_range=[0, T], y_range=[0, 1500000])
plot_1.line('x', 'y', source=source_1, line_width=3, line_alpha=0.6)
tab1 = Panel(child=plot_1, title="Recovered")

plot_2 = figure(plot_height=400, plot_width=400, title="SEIR model",
                tools="crosshair,pan,reset,save,wheel_zoom",
                x_range=[0, T], y_range=[0, 20000])
plot_2.line('x', 'y', source=source_2, line_width=3, line_alpha=0.6)
tab2 = Panel(child=plot_2, title="Infected")

tabs = Tabs(tabs=[tab1, tab2])


def update_data(attrname, old, new):

    # Get the current slider values
    incub = ip.value
    sel_val = select.value
    active = button.active
    # Generate the new curve
    x = np.linspace(0, T, 100000)

    y = model_solution(x,
                       r_0=r0_list[active],
                       comp='R',
                       population=pop_dict[sel_val],
                       incubation_period=incub)
    source_1.data = dict(x=x, y=y)
    y = model_solution(x,
                       r_0=r0_list[active],
                       comp='I',
                       population=pop_dict[sel_val],
                       incubation_period=incub)
    source_2.data = dict(x=x, y=y)


for w in [select, ip]:
    w.on_change('value', update_data)
button.on_change('active', update_data)

# Set up layouts and add to document
inputs = column(ip, button, select)

curdoc().add_root(row(inputs, tabs, width=800))
curdoc().title = "SEIR"
