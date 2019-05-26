"""

    Author: Udi Yosovzon
    Email: udiyosovzon@gmail.com
    ID: 308063437
    
    Description: This python code evaluate a specific inegral using Trapez method and Richardson extrapolation method.
    			 It then compares to numpy result
                 
    Execution: In cmd, navigate to this file directory, then type: "python Q1.py"

"""


import numpy as np


def f(x):
    
    """ defining the desired function
    """
    return np.exp(-(x**2))


def trapez(f, a, b, n=20):
    
    """ f - function to integrate
        a - lower limit
        b - upper limit
        n - number of segments
    """
    h = (b-a)/n    # step
    internal_points = [2*f(a + i*h) for i in range(1,n)]
    total_sum = (h/2)*(f(a) + sum(internal_points) + f(b))
    
    return total_sum


def print_solution(sol, sol_name):
    
    """
    """
    print(f"\n{sol_name} solution: {sol}")
    print(f"Difference from NumPy solution: {np.abs(sol - np_sol)}")
    print("-"*100)



a = 0    # lower limit
b = 2    # upper limit
x = np.arange(0,2,0.000001)
np_sol = np.trapz(f(x),x)
print(f"\nAnalytic solution using NumPy: {np_sol}")
print("-"*100)

trapez_sol = trapez(f,a,b)
print_solution(trapez_sol, "Trapez")

richardson_sol = (4*trapez(f,a,b) - trapez(f,a,b,10))/3
print_solution(richardson_sol, "Richardson")