"""

	Author: Udi Yosovzon
	
		
	Description: Helper function to calculate cubic spline
	Execution: No need to execute directly

"""

import numpy as np

def solve_banded_system(x,y,h):
    
    """ takes in lists of x,y values and h values and cosntruct a banded system (tridiagonal matrix). 
        solves the system and returns result
    """
    n_unknowns = len(x) - 2    # number of unknowns
    
    A = np.zeros((n_unknowns,n_unknowns))
    v = np.zeros(n_unknowns)
        
    for i in range(n_unknowns):
        A[i][i] = 2*(h[i] + h[i+1])
        v[i] = 6*((y[i+2]-y[i+1])/h[i+1] - (y[i+1]-y[i])/h[i])
    
    for i in range(n_unknowns-1):
        A[i][i+1] = h[i+1]
        A[i+1][i] = h[i+1]
    
    return np.linalg.solve(A, v)

def construct_splines(fx_mapping):
    
    """ takes in a sequence of (x,y) points and for each section between two points calculates a cubic spline
    """
    n_points = len(fx_mapping)    # number of points
    n_splines = n_points - 1    # number of splines
    
    x = [item[0] for item in fx_mapping]    # get x values
    y = [item[1] for item in fx_mapping]    # get y values
    h = [x[i+1] - x[i] for i in range(n_splines)]    # calculate h values from x values
    
    # for each point calculate a value (store in m) - this value will be later used to calculate all coefficients
    m = np.zeros(n_points)
    m[1:-1] = solve_banded_system(x,y,h)
    
    # create a structured array to store spline coefficients and section lim[i]ts
    splines = np.zeros(n_splines, dtype=[("start", np.float64), ("end", np.float64), ("coefs",np.float64, 4)])
    
    # calculate coefficients using predefined, well known formulas
    for i in range(n_splines):
        
        c1 = (m[i+1]-m[i])/(6*h[i])
        c2 = (x[i+1]*m[i] - x[i]*m[i+1])/(2*h[i])
        c3 = (x[i]**2*m[i+1] - x[i+1]**2*m[i])/(2*h[i]) + (y[i+1]-y[i])/h[i] + (m[i]-m[i+1])*h[i]/6
        c4 = (x[i+1]**3*m[i]-x[i]**3*m[i+1])/(6*h[i]) + (x[i+1]*y[i]-x[i]*y[i+1])/h[i] + (x[i]*m[i+1]-x[i+1]*m[i])*h[i]/6
        
        splines[i]["start"] = x[i]
        splines[i]["end"] = x[i+1]
        splines[i]["coefs"] = [c1,c2,c3,c4]
        
    return splines

def make_func(coefs):
    
    """ coefs is a list/array of a polynomial coefficients.
        returns a convenient function that calculates function value at given point x
    """
    n = len(coefs)
    
    def f(x):
        return sum([coefs[i]*x**(n-1-i)  for i in range(n)])
    
    return f

def calculate_cubic_spline(fx_mapping, h):
    
    """
    """
    splines = construct_splines(fx_mapping)
    
    x = []
    f = [] # unify all section splines into one big list

    for spl in splines:
        start = spl["start"]
        end = spl["end"]
        
        if spl==splines[-1]:    # it this is the last section, include last point
            x_spl = np.arange(start,end+1,h)
        else:
            x_spl = np.arange(start,end,h)
        
        f_spl = make_func(spl["coefs"])
        
        x += list(x_spl)
        f += list(f_spl(x_spl))
        
    n = len(f)
    return {x[i]:f[i] for i in range(n)}