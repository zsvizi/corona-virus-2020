import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as c
import matplotlib.ticker as ticker
from decimal import Decimal

from risk import get_heatmap

'''
Hello
heatmap.py plots the output of risk.py as a heatmap graph with the given parameters.

The risk is a function of 3 variables: C, theta and r_loc. We follow this order in this script.
        - i) fix one variable: param = (YOUR SELECTED VARIABLE) 
                    ex. param = 0, that is "C"  (as Python idexes from 0) 
                    or param = 1, that is "theta"
                    or param = 2, that is "r_loc"
        - ii) say, you selected param = 2 (r_loc) to be fixed, now you can plot as many heatmaps as many       r_loc values you provide to the program in r_loc = [(YOUR LIST OF FIXED PARAMETERS)]
                    ex. r_loc = np.array([1.3,1.8,2.2,2.6])
        - iii) specify the range of the other two variables you want to examine
                    ex. c = np.arange(100000, 1100000, 10000) 
                    and theta = np.arange(0, 0.000255, 0.000001)
        - iv) The program saves the heatmaps in a folder named after your selected variable in the             "Heatmaps" folder
'''
#----------------------------------------------------        


# Give the values here:
param = 2           #param = 0 ( = C); param = 1 ( = theta); param = 2 ( = r_loc)

c = np.arange(100000, 1100000, 10000)
theta = np.arange(0, 0.000255, 0.000001)
r_loc = np.array([1.05, 1.3, 1.8, 2.2, 2.6])
#r_loc = np.arange(1, 2.5, 0.1)

#----------------------------------------------------        

# Don't touch this :p

coords = [0,1,2]
coords.pop(param)
vals = [c, theta, r_loc][param]
name_of_coords = ["C", "$\\theta$", "$R_{loc}$"]

def get_data(c, theta, r_loc):
    heatmap, _ = get_heatmap(c, theta, r_loc)
    c = heatmap["r_stars"]
    theta = heatmap["theta"]
    r_loc = heatmap["r_locs"]
    coord_list = heatmap["heatmap"]
    
    return [[c, theta, r_loc], coord_list]


def preprocess():
    variables, raw_data = get_data(c, theta, r_loc)
    data = []
    
    for value in vals:
        temp_data = []
        for element in raw_data:    # goes through raw_data
            # selects required values
            if abs(element[param] - value) < 10 ** (-15):  
                # ads it to a list
                temp_data.append([element[coords[0]], element[coords[1]], element[-1]])
        data.append(temp_data)
    data = np.array(data) 
    return [variables, data]


def plot_heatmap(variables, vals, dat,i,foldr):
    
    X = variables[coords[0]]
    Y = variables[coords[1]]
    
    dat = np.array(dat[:, -1]).reshape(len(X), len(Y))
    Z = dat.T

    fig, ax = plt.subplots()
    plt.contourf(X, Y, Z, 50, cmap='Reds', alpha=1)

    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())   
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    ax.set_title(name_of_coords[param]+' = ' + str(vals[i]), fontsize=10)
    plt.xlabel(name_of_coords[coords[0]])
    plt.ylabel(name_of_coords[coords[1]])

    plt.clim(0,1)
    plt.colorbar()

    plt.tight_layout()
    plt.savefig(foldr+"/"+str(vals[i])+".png", dpi=300, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None)
    
    plt.clf()

#----------------------------------------------------        
    

def makedir():
    names = ["C", "Theta", "R_loc"]
    foldr = "./Heatmaps/" + names[param]
    try:
        os.makedirs(foldr)
    except:
        foldr = foldr+"-"+datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        os.makedirs(foldr)
    return foldr


#----------------------------------------------------      


def main():
    foldr = makedir()       
    variables, data = preprocess()

    for i in range(len(data)):
        plot_heatmap(variables, vals, data[i], i, foldr)
    
    

if __name__ == "__main__":
    main()

