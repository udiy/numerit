"""

	Author: Udi Yosovzon
	Email: udiyosovzon@gmail.com
	ID: 308063437
	
	Description: This python code solves a system of linear eqations with three methods:
					1) Gauss elimanation
					2) LU decomposition
					3) Gauss Seidel
				 The user is presented a menu from which they can choose one the three above
	Execution: In cmd, navigate to this file directory, then type: "python Q1.py"

"""


import numpy as np

A = np.array([
        [3,-3,2,-4],
        [-2,-1,3,-1],
        [5,-2,-3,2],
        [-2,4,1,2]
    ], dtype=np.float64)

b = np.array([
    [7.9],
    [-12.5],
    [18],
    [-8.1]
], dtype=np.float64)


def division_or_warning(a,b):
    
    """ divide a by b (a/b) and return the result
        if b is zero, do nothing, and print a warning
    """
    if b != 0:
        return a/b
    else:
        print("Division by zero is not possible!")
        return None


def is_triangular(A):
    
    """ check if matrix A is triangular than return if it is upper or lower triangular
    """
    m,n = A.shape
    
    # verify number of rows and number of columns are the same
    if m != n:
        print("Your matrix is not a square matrix! Aborting substitution!")
        return
    else:
        lower = True
        upper = True
        # iterate n-1 iterations and count number of zeros
        for i in range(0,n-1):
            if lower:
                # make sure the first i elements of the i+1 row are zero
                num_of_zeros_lower = np.count_nonzero(A[i+1][:i+1])
                if num_of_zeros_lower != 0:
                    lower = False
            if upper:
                # or the last n-(i+1) elements of the i row are zero
                num_of_zeros_upper = np.count_nonzero(A[i][i+1:])
                if num_of_zeros_upper != 0:
                    upper = False
            # if we found out the matrix is not upper or lower triangular, break out from the for loop
            if not (upper or lower):
                break
        
        if lower and not upper:
            return "lower"
        elif upper and not lower:
            return "upper"
        elif lower and upper:
            print("Your matrix is diagonal")
            return None
        print("Could not classify your matrix as either lower or upper triangular")
        return None
        
        
def substitution(A,b):
    
    """ perform forward or backward substitution for a triangular matrix.
        Upper triangular - all elements in the upper triangle are 0 - forward substitution
        Lower triangular - all elements in the lower triangle are 0 - backward substitution
    """
    # x is variables vector, has the same dimensions as b
    x = np.ones_like(b)
    m = x.shape[0]
    
    # verify this is triangular matrix then determine subs type
    tr_type = is_triangular(A)
    
    # determine variables start, end and step for the for loop to calculate x vector
    # if this is an upper triangular matrix set the the variables for forward substitution
    if tr_type == "upper":
        start = 0
        end = m
        step = 1
    # if this is an lower triangular matrix set the the variables for backward substitution
    elif tr_type == "lower":
        start = m - 1
        end = -1
        step = -1
    else:
        print("Cannot perform substitution with your matrix")
        return
    
    # calculate xi with the substitution formula
    for i in range(start, end, step):
        _sum = sum([A[i][j] * x[j] for j in range(start, i, step)])
        
        normalize_factor = division_or_warning(1,A[i][i])        
        if normalize_factor:
            x[i] = (normalize_factor * (b[i] - _sum))
            
        
    return x


def gauss_elimanation(A,b):
    
    """ A - coefficient matrix
        b - solution vector
    """
    # add solution vector b as another column at the end of matrix A
    C = np.c_[A, b]
    
    # get number of rows in C
    m = C.shape[0]
    
    for i in range(m):
        pivot_element = C[i][i]
        
        # for the rows below that row, subtract this row, so the lower triangle of the matrix is all zeros
        for j in range(i+1, m):
            normalize_factor = division_or_warning(C[j][i],pivot_element)
            if normalize_factor:
                C[j] = C[j] - (normalize_factor)*C[i]
                
    lower_tr_A = C[:,:-1]
    new_b = C[:,-1:]
    return lower_tr_A, new_b


def lu_decomposition(A):
    
    """ decompose matrix A to L and U
    """
    # get number of rows in A, assuming a is n by n matrix
    n = A.shape[0]
    
    # initialize L as n by n I matrix
    L = np.eye(n)
    
    U = A.copy()
    
    for i in range(1,n):
        
        # for row index i in A (zero indexed) we need i zeros to get U, and we need to calculate i elements for L
        # j is the index of the column and U[j][j] are the elements on the diagonal, for the j-th column
        for j in range(i):
            uij = division_or_warning(U[i][j],U[j][j])
            if uij:
                L[i][j] = uij
                U[i] = U[i] - uij*U[j]
            
    
    return L,U


def pivot_matrix_rows(A, b):
    
    """ A is a square matrix
    """
    m,n = A.shape
    
    A_pivoted = np.zeros_like(A)
    b_pivoted = np.zeros_like(b)
    
    # for each row in matrix A find the j index (which column) of its max element
    for i in range(m):
        max_element_index = np.abs(A[i]).argmax()
        
        # make sure this index of max element is not in A pivoted already
        if np.count_nonzero(A_pivoted[max_element_index]) == 0:
            # inset row i from A to A_pivoted so the diagonal element is the max element
            A_pivoted[max_element_index] = A[i]
            b_pivoted[max_element_index] = b[i]
        else:
            print("This row is already taken!")
    
    return A_pivoted, b_pivoted


def gauss_seidel(A, b, x_init, epsilon=10**-9, i_num=1):
    
    """ x_init is the initial guess for x
    """
    # pivot A rows
    A_piv, b_piv = pivot_matrix_rows(A, b)
    
    m = x_init.shape[0]
    x = x_init.copy()
    
    for i in range(m):
        
        _sum = sum([(A_piv[i][j] * x[j]) for j in range(m) if j!=i])
        
        x[i] = (b_piv[i] - _sum) / (A_piv[i][i])
    
    # check convergence
    max_diff = np.abs(x - x_init).max()
    if max_diff < epsilon:
        return x
    elif i_num > 1000:
        print(f"Reached maximal number of iterations: {i_num}")
        return
    else:
        return gauss_seidel(A, b, x, epsilon, i_num=i_num+1)




methods = {
        1: "Gauss elimanation",
        2: "LU decomposition",
        3: "Guass Seidel"
    }

def print_menu():
    
    print("Choose from the following, or type '#' to exit:\n")
    for i, method in methods.items(): 
        print(f"\tType {i} for {method}")
    return input()


if __file__ == "Q1.py":
	while True:
	    x = 0
	    met = print_menu()

	    if met == "#":
	        break

	    met = int(met)
	    if met in methods.keys():
	        print(f"You chose {methods[met]}")
	        
	        # Solving with gauss elimanation
	        if met == 1:
	            A_elimanated, b_elimanated = gauss_elimanation(A,b)
	            x = substitution(A_elimanated, b_elimanated)
	            
	        # Solving with lu decomposition
	        elif met == 2:
	            L, U = lu_decomposition(A)
	            z = substitution(L, b)
	            x = substitution(U, z)
	            
	        # Solving with Gauss Seidel
	        elif met == 3:
	            x_init = np.array([[0], [0], [0], [0]], dtype=np.float64)
	            x = gauss_seidel(A, b, x_init)
	        
	        print(f"Solution found:\n\n{x}\n\n")
	    
	    else:
	        print(f"{met} is not an option! Try again")