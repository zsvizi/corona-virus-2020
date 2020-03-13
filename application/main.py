# To run this code, create a Bokeh server via
# `bokeh serve --show application`
# from folder ./corona-virus-2020

import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, RadioButtonGroup, Tabs, Panel, Select
from bokeh.plotting import figure

from .utils import EpidemicModel, solve_model


# Set up global variables
hungary_population = 9667595
t_range = 200
pop_dict = {"Magyarország": 1,
            "Budapest": 5,
            "Debrecen": 100}
r0_list = [2.1, 2.8]
y_range = 200


def model_solution(t, r_0, comp, city):
    # Set time intervals
    infectious_period = 3.3
    incubation_period = 5.1
    exposed_init = np.array([2.0, 0.0])/city
    infected_init = np.array([0.0, 0.0, 0.0])/city
    t_star = 50
    city_population = hungary_population/city

    beta = r_0 / (city_population * infectious_period)
    alpha = 1 / incubation_period
    gamma = 1 / infectious_period
    params = np.array((beta, alpha, gamma))
    # Set initial values
    init_values = {"e0": exposed_init,
                   "i0": infected_init,
                   "r0": [0],
                   "c0": [0]
                   }
    not_susceptible = sum([item for sublist in init_values.values() for item in sublist])
    init_values.update({"s0": [city_population - not_susceptible]})

    model = EpidemicModel(init_values, t_star=t_star, r_0=r_0)
    x0 = np.array(model.get_initial_values())
    # Solve epidemic model
    solution = solve_model(t, x0, params, model)
    sol = []
    if comp == 'C':
        sol = solution[:, -1]
    elif comp == 'I':
        sol = solution[:, -3] + solution[:, -4] + solution[:, -5]
    return sol


# Set up widgets
# ip = Slider(title="Incubation period", value=5.1, start=5.0, end=14.0, step=0.1)
button = RadioButtonGroup(labels=["R0=2.1", "R0=2.8"], active=0)
select = Select(title="Város:", value="Magyarország", options=["Magyarország", "Budapest", "Debrecen"])

# Initial plot
x = np.linspace(0, t_range, 100000)
y = model_solution(x,
                   comp='C',
                   r_0=r0_list[button.active],
                   city=pop_dict[select.value])
source_1 = ColumnDataSource(data=dict(x=x, y=y))
y = model_solution(x,
                   comp='I',
                   r_0=r0_list[button.active],
                   city=pop_dict[select.value])
source_2 = ColumnDataSource(data=dict(x=x, y=y))


# Set up plot
plot_1 = figure(plot_height=600, plot_width=600, title="Járványterjedési modell ",
                tools="crosshair,pan,reset,save,wheel_zoom",
                x_range=[0, t_range], y_range=[0, y_range],
                x_axis_label="Március 1. óta eltelt napok száma", y_axis_label="Esetszám")
plot_1.xaxis.axis_label_text_font_size = "15pt"
plot_1.yaxis.axis_label_text_font_size = "15pt"

plot_1.line('x', 'y', source=source_1, line_width=3, line_alpha=0.6)
tab1 = Panel(child=plot_1, title="Kumulatív fertőzött")

plot_2 = figure(plot_height=600, plot_width=600, title="Járványterjedési modell",
                tools="crosshair,pan,reset,save,wheel_zoom",
                x_range=[0, t_range], y_range=[0, y_range/7],
                x_axis_label="Március 1. óta eltelt napok száma", y_axis_label="Esetszám")
plot_2.line('x', 'y', source=source_2, line_width=3, line_alpha=0.6)
plot_2.xaxis.axis_label_text_font_size = "20pt"
plot_2.yaxis.axis_label_text_font_size = "20pt"
tab2 = Panel(child=plot_2, title="Fertőzöttek")

tabs = Tabs(tabs=[tab1, tab2])


def update_data(attrname, old, new):

    # Get the current slider values
    sel_val = select.value
    active = button.active
    # Generate the new curve
    x = np.linspace(0, t_range, 100000)

    y = model_solution(x,
                       r_0=r0_list[active],
                       comp='C',
                       city=pop_dict[sel_val])
    source_1.data = dict(x=x, y=y)
    y = model_solution(x,
                       r_0=r0_list[active],
                       comp='I',
                       city=pop_dict[sel_val])
    source_2.data = dict(x=x, y=y)


for w in [select]:
    w.on_change('value', update_data)
button.on_change('active', update_data)

# Set up layouts and add to document
inputs = column(button, select)

curdoc().add_root(row(inputs, tabs, width=1200))
curdoc().title = "SEIR"
