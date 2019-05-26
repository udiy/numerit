"""

	Author: Udi Yosovzon
	
		
	Description: This python code uses two methods of interpolation to evaluate a specific function.
				 First method is "Direct (Yeshara) method", it involves solving a system of liner eqations to find coefficients for function to evaluate x.
				 Second method is Lagrange method, it calculates weights for every data point and calculate value at x.
				 
	Execution: In cmd, navigate to this file directory, then type: "python Q1.py"

"""


import numpy as np


# each key value pair is x as the key, and fx as the value
fx_mapping = { -1:3, 0:1, 2:2 }
x = -0.5



def construct_vandermonde_matrix(x_values):
    
    """ create a vandermonde matrix
    """
    n = len(x_values)
    x_values = list(x_values)
    
    vandermonde = np.zeros((n,n))
    for i in range(n):
        x = x_values[i]
        for j in range(n):
            vandermonde[i][j] = x**j
    
    return vandermonde



def yeshara(fx_mapping, x):
    
    """ fx_mapping - dict of x:fx values
        x - point to evaluate function at
    """
    x_values = fx_mapping.keys()
    fx_values = np.array(list(fx_mapping.values()))
    
    # calculate coefficient for polynom
    vandermonde_matrix = construct_vandermonde_matrix(x_values)
    coef_vec = np.linalg.solve(vandermonde_matrix, fx_values)
    
    # evaluate polynom at x
    n = coef_vec.shape[0]
    poly_at_x = sum([coef_vec[i]*(x**i) for i in range(n)])
    
    return poly_at_x



def calculate_lagrange_weight(x_values, i, x):
    
    """ x_values - list of x_values
        i - index of L
        x - point of evaluation of L 
    """
    n = len(x_values)
    L = 1
    
    for j in range(n):
        if j != i:
            xj = x_values[j]
            xi = x_values[i]
            L *= (x - xj)/(xi - xj)
            
    return L



def lagrange(fx_mapping, x):
    
    """ fx_mapping - dict of x:fx values
        x - point to evaluate function at
    """
    x_values = list(fx_mapping.keys())
    fx_values = list(fx_mapping.values())
    
    n = len(fx_mapping)
    poly_at_x = []
    
    for i in range(n):
        
        # calculate L - lagrange weight at i, x
        L = calculate_lagrange_weight(x_values, i, x)
        
        # multiply by fx and append to list
        poly_at_x.append(L*fx_values[i])
        
    return sum(poly_at_x)


if __file__ == "Q1.py":
	result_yeshara = yeshara(fx_mapping, x)
	result_lagrange = lagrange(fx_mapping, x)
	print(f"""
	Given interpolation points {list(fx_mapping.items())},
	we want to construct a 2nd order polynomial f
	and calculate f(-0.5) in two different methods.
	Results:
	            Yeshara: f({x}) = {result_yeshara}
	            Lagrange: f({x}) = {result_lagrange}
	""")