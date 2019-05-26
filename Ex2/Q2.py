"""

	Author: Udi Yosovzon
	
		
	Description: This python code computes the solution to a specific system of eqations
				 It then prints the solution to screen.
	Execution: In cmd, navigate to this file directory, then type: "python Q2.py"

"""

import numpy as np


eq1 = np.array([
    [0,0,-52], # x to the power of 1
    [4,4,-19] # x to the power of 0
])

eq2 = np.array([
    [0, 0, 169], # x to the power of 2
    [0, 0, -111], # x to the power of 1
    [3, -10, 0] # x to the power of 0
])

sys_eq = [
    eq1,
    eq2
]


starting_point = (-0.01, -0.01)


def evaluate_polynomial_2d(poly2d, point):
    
    """ poly - polynomial with 2 variables, will be represented by a two dimensional coefficient matrix.
        row index is the order of x, column index is the order of y.
        For example: 
            [
                [0,0,-52], # x to the power of 1
                [4,4,-19] # x to the power of 0
            ]
        point - a tuple (x,y). For example: (-0.01, -0.01)
    """
    x, y = point
    poly2d = np.array(poly2d)
    m, n = poly2d.shape
    _sum = 0
    
    for i in range(m):
        x_order = m-(i+1)
        for j in range(n):
            y_order = n-(j+1)
            coef = poly2d[i][j]
            _sum += coef*(x**x_order)*(y**y_order)
    
    return _sum
           



def evaluate_partial_derivatives(poly2d, point):
    
    """
    """
    x, y = point
    x1 = x - 0.0001
    y1 = y - 0.0001
    
    # polynomial value at x,y
    pxy = evaluate_polynomial_2d(poly2d, point)
    # polynomial value at x1,y
    px1y = evaluate_polynomial_2d(poly2d, (x1,y))
    # polynomial value at x,y1
    pxy1 = evaluate_polynomial_2d(poly2d, (x,y1))
    
    # partial derivative by x, while y remains is constant
    x_derivative = (pxy - px1y)/(x - x1)
    # partial derivative by y, while x remains is constant
    y_derivative = (pxy - pxy1)/(y - y1)
    
    return x_derivative, y_derivative 



def newton_raphson_2d(equations, point, epsilon, i=1):
    
    """
    """
    eq_value_at_point = []
    
    # evaluate system of eqations at point
    for eq in equations:
        # evaluate every eqation (polynomial) at point and append result to eq_value_at_point
        eq_value = evaluate_polynomial_2d(eq, point)
        eq_value_at_point.append(eq_value)
    
    # if evalutation is close enough to zero (less than epsilon) return point
    if all([np.abs(val) <= epsilon for val in eq_value_at_point]):
        print(f"\tnumber of iterations: {i}")
        return point
    # else use newton formula to calucalte next point
    else:
        eq_derivative_at_point = []
        # evaluate partial derivaites of every eqation
        for eq in equations:
            eq_derivaitves = evaluate_partial_derivatives(eq, point)
            eq_derivative_at_point.append(eq_derivaitves)
        
        f, g = eq_value_at_point
        fx, fy = eq_derivative_at_point[0]
        gx, gy = eq_derivative_at_point[1]
        
        denominator = (fx*gy - gx*fy)
        
        if denominator != 0:
            x,y = point
            x_new = x - (f*gy - g*fy)/denominator
            y_new = y - (g*fx - f*gx)/denominator
            # call newton_raphson_2d again with new point
            return newton_raphson_2d(equations, (x_new, y_new), epsilon, i=i+1)
        else:
            print(f"\tnumber of iterations: {i}")
            return point



solution = newton_raphson_2d(sys_eq, starting_point, epsilon=10**-4)
print(f"\n\tThe solution found with newton_raphson_2d:\n\n\t{solution}")
print("-"*60)