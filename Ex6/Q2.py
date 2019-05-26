"""

    Author: Udi Yosovzon
    Email: udiyosovzon@gmail.com
    ID: 308063437
    
    Description: This python code evaluates a specific integral using the following methods:
    			1) Trapez
    			2) Simpson
    			3) Romberg
    			4) Quadrature
    			It then compares to the result calculated with numpy
                 
    Execution: In cmd, navigate to this file directory, then type: "python Q2.py".

"""

import re
import numpy as np

def f(x):
    
    """ defining the desired function
    """
    return x*np.exp(2*x)

def trapez(f, a, b, n=20):
    
    """ evaluate integral using trapez method
        f - function to integrate
        a - lower limit
        b - upper limit
        n - number of segments
    """
    h = (b-a)/n    # step
    internal_points = [2*f(a + i*h) for i in range(1,n)]
    total_sum = (h/2)*(f(a) + sum(internal_points) + f(b))
    
    return total_sum

def simpson(f, a, b, n=20):
    
    """ evaluate integral using Simpson 1/3 method
        f - function to integrate
        a - lower limit
        b - upper limit
        n - number of segments
    """
    h = (b-a)/n    # step
    _sum = 0
    
    for i in range(n):
        start = a + i*h    # starting point of the ith segment
        end = start + h    # ending point of the ith segment
        middle = (start+end)/2    # middle point of the ith segment
        _sum += (f(start) + 4*f(middle) + f(end))
        
    return (h/6)*_sum


# Romberg functions

def while_generator(h):
    
    """ an helper function to use for a list comprehension
        it returns a list generator every item a list in the list is number of segments to use for evaluation of
        simpsons/trapez method. For example: if b-a=4 than it returns [1, 2, 4, 8, 16]
    """
    n = 1    # number of segments
    while True:
        yield n
        n = n*2    # increase number of semgments by a factor of 2
        if n > h**2:    # stop if number of segments if greater than h**2 (happens when h=1/(b-a))
            break

def romberg(f, a, b, method):
    
    """ evaluate integral using Romberg method
        f - function to integrate
        a - lower limit
        b - upper limit
        method - either simpson or trapez
    """
    h = b - a
    I = [method(f,a,b,i) for i in while_generator(h)]    # store method evaluations in this list
    k = 1    # index for Romberg iteration
    
    while True:
        print(I)
        n = len(I)
        if n == 1:    # if we got to the final evaluation break the while loop and return the last evaluation
            return I[0]
        
        I = [((4**k)*I[i+1] - I[i])/(4**k - 1) for i in range(n-1)]    # use Richardson formula to evaluate next level
        k += 1


# Quad functions

def make_point_tuple(point, negative=False):
    
    """ takes in a list of strings, first item is point value, second is weight
        returns a tuple of floats
    """
    point_value = float(point[0])
    weight = float(point[1])
    if negative:
        return (-point_value, weight)
    return (point_value, weight)

def point_string(point):
    
    """ takes in a point, weight represantation as string and returns a list of one tuple 
        or two tuples if there a is a plus-minus sign
    """
    point = point.split()
    if point[0] == "±":            
        return [make_point_tuple(point[1:]), make_point_tuple(point[1:], negative=True)]
    else:
        return [make_point_tuple(point)]

def points_string(points):
    
    """ takes in a string representing points and weights for each point, returns a list of tuples
    """
    points = points.strip()    # trim spaces on the edges
    points = points.split("\n")    # split points to a list
    points_tuples = []
    for point in points:
        points_tuples += point_string(point)
        
    return points_tuples

quad_values = """
2 ± 0.57735026 1.0
3 0.0 0.88888889
± 0.77459667 0.55555555
4 ± 0.33998104 0.65214515
± 0.86113631 0.34785485
5 0.0 0.56888889
± 0.53846931 0.47862867
± 0.90617985 0.23692689
6 ± 0.23861918 0.46791393
± 0.66120939 0.36076157
± 0.93246951 0.17132449
7 0.0 0.41795918
± 0.40584515 0.38183005
± 0.74153119 0.27970539
± 0.94910791 0.12948497
8 ± 0.18343464 0.36268378
± 0.52553241 0.31370665
± 0.79666648 0.22238103
± 0.96028986 0.10122854
10 ± 0.14887434 0.29552422
± 0.43339539 0.26926672
± 0.67940957 0.21908636
± 0.86506337 0.14945135
± 0.97390653 0.06667134
"""

def quadrature(f, a, b, quad_values):
    
    """
    """
    # handling quad_values string, make it a useable data structure
    idx = re.findall("\\n\d{1,2}", quad_values)
    points = re.split("\\n\d{1,2}",quad_values)[1:]
    n_cases = len(idx)
    idx = [int(idx[i][1:]) for i in range(n_cases)]
    quad_values_dict = {idx[i]: points_string(points[i]) for i in range(n_cases)}
    
    for n,pts in quad_values_dict.items():
        I = sum([((b-a)/2)*p[1]*f(((b-a)/2)*p[0] + (b+a)/2) for p in pts])
        print(f"Quad {n} points: {I}")
        
    return I

# Main

def print_solution(sol, sol_name):
    
    """
    """
    print(f"\n{sol_name} solution: {sol}")
    print(f"Difference from NumPy solution: {np.abs(sol - np_sol)}")
    print("-"*100)


a = 0    # lower limit
b = 4    # upper limit
x = np.arange(a,b,0.000001)
np_sol = np.trapz(f(x),x)
print(f"\nAnalytic solution using NumPy: {np_sol}")
print("-"*100)

trapez_sol = trapez(f,a,b)
print_solution(trapez_sol, "Trapez")

simpson_sol = simpson(f,a,b)
print_solution(simpson_sol, "Simpson")

print("\nRomberg")
romberg_sol = romberg(f,a,b,simpson)
print_solution(romberg_sol, "Romberg")

print("\nQuad")
quadrature_sol = quadrature(f,a,b,quad_values)
print_solution(quadrature_sol, "Quad")