"""

    Author: Udi Yosovzon
    Email: udiyosovzon@gmail.com
    ID: 308063437
    
    Description: This python code creates a Lagrange polynomial for the function cos(4*pi*x), than compare the results
                 of the interpolation and the function.
                 
    Execution: In cmd, navigate to this file directory, then type: "python Q3.py". IMPORTANT: this file has to be in the same folder as Q1

"""


import numpy as np
import matplotlib.pyplot as plt
import sympy
from sympy.abc import pi as pi_sym
from Q1 import lagrange



def f(x):

    """ cos(4*pi*x) function
    """
    pi = np.pi
    cos = np.cos(4*pi*x)
    round_result = np.around(cos, decimals=5)
    return round_result


# sample 4 points from function f, evenly spaced in the range [0,0.5]
x_values = np.linspace(0,0.5,4)
fx_mapping = {x:f(x) for x in x_values}
fx_values = fx_mapping.values()


# ask user to input a number and calculate and print results
while True:
    
    input_x = input("Enter a number in the range of [0,0.5] or type '#' to quit: ")
    if input_x == "#":
        break
    else:
        input_x = float(input_x)
        print(f"You entered: {input_x}")

        res_interpolation = lagrange(fx_mapping,0.5)
        res_function = f(input_x)

        print(f"""
        Results\tInterpolation\tFunction
               \t{res_interpolation}          \t{res_function}
        """)



# optionally plot functions
input_plot = input("Do you want to plot functions? (y/n) ")
if input_plot.lower() == "y":
    x = np.linspace(0,0.5,num=1000)
    plt.figure(figsize=(16,9))
    plt.plot(x,f(x),label=f"numpy cos(4{sympy.pretty(pi_sym)}x) function")
    plt.plot(x,lagrange(fx_mapping,x),label="Lagrange interpolation")
    plt.plot(x_values,fx_values,"o",c="red", markersize=15)
    plt.axhline(linewidth=1, color="lightgrey")
    plt.legend(prop={'size': 16})
    plt.show()