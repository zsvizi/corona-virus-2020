import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as c
import matplotlib.ticker as ticker

from source.risk import get_heatmap

# DIMENSIONS: r_star (final size at t_star, later maybe t_star) x connectivity (theta) x r_loc (~z)
param = 0
coords = [1, 2]
vals = [1.05]
filename = "..\\data\\teszt_2.txt"


def get_data(path):
    max_number_summands = 100
    read = False
    if read:
        coord_list = read_data(path)
        r_loc = coord_list[0]
        theta = coord_list[1]
        t_star = coord_list[2]
        del coord_list[:3]
    else:
        coord_list, t_star, theta, r_loc = get_heatmap(max_number_summands)
    
    return [[r_loc, theta, t_star], coord_list]


def read_data(path):
    # beolvassa az aktuális (nr.txt) fájlt, soronként feldarabolja és a ccord_list-be teszi
    with open(path, "r") as coord_file:
        coord_list = coord_file.read().splitlines()
    # végigmegy a coord_list elemein és a "," mentén feldarabolja
    for number in range(len(coord_list)):
        row = coord_list[number].split(",")
        coord_list[number] = []
        for elem in row:
            coord_list[number].append(float(elem))
    return coord_list


def preprocess(fname):
    vars, raw_data = get_data(fname)
    data = []
    for value in vals:
        temp_data = []
        for element in raw_data:    # végig megy a teljes nyers adatsoron
            if abs(element[param] - value) < 10 ** (-4):  # kiválogatja a megadott paraméter szerint a kért értékeket
                # hozzáadja egy listához
                temp_data.append([element[coords[0]], element[coords[1]], element[-1]])
        data.append(temp_data)
    data = np.array(data)
    return [vars, data]


def plot_heatmap(dat,i,foldr):
    
    X = variables[coords[0]]
    Y = variables[coords[1]]
    """
    for value in Y:
        lst = []
        for row in dat:
            if row[1] == value:
                lst.append(row[-1])
        Z.append(lst)
    """
    dat = np.array(dat[:, -1]).reshape(len(X), len(Y))
    Z = dat.T

    fig, ax = plt.subplots()
    contours = plt.contour(X, Y, Z, [0.3,0.4,0.5], colors='#000000')
    plt.clabel(contours, inline=True, fontsize=8)
    plt.contourf(X, Y, Z, 10, cmap='RdGy', alpha=0.5)

    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())   
    
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    
    if coords[0] == 0: X_name = "local reproduction number"
    elif coords[0] == 1: X_name = "Theta"
    elif coords[0] == 2: X_name = "R*"
    
    if coords[1] == 0: Y_name = "local reproduction number"
    elif coords[1] == 1: Y_name = "Theta"
    elif coords[1] == 2: Y_name = "R*"

    plt.xlabel(X_name)
    plt.ylabel(Y_name)
    
    plt.colorbar()
    plt.savefig(foldr+"/"+str(vals[i])+".png", dpi=300, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None)
    plt.clf()


def makedir():
    if param == 0: foldr = "./R_loc"
    elif param == 1: foldr = "./theta"
    elif param == 2: foldr = "./t_star"
    else: foldr = "valami_baly_van!4444!!"+str(vals)
    
    try:
        os.makedirs(foldr)
    except:
        dt = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
        foldr = foldr+"-"+dt
        os.makedirs(foldr)
    
    return foldr

foldr = makedir()

variables, data = preprocess(filename)
for i in range(len(data)):
    plot_heatmap(data[i],i,foldr)
