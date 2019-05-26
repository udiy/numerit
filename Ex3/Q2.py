"""

	Author: Udi Yosovzon
	Email: udiyosovzon@gmail.com
	ID: 308063437
	
	Description: This python code finds the inverse of matrix A using the method LU decomposition
	
	Execution: In cmd, navigate to this file directory, then type: "python Q2.py". 
               IMPORTANT - this file must be in the same directory as Q1.py because it has code dependencies!

"""


import numpy as np
from Q1 import lu_decomposition, substitution, gauss_elimanation

A = np.array([
    [4,8,4,0],
    [1,4,7,2],
    [1,5,4,-3],
    [1,3,0,-2]
], dtype=np.float64)

L, U = lu_decomposition(A)

print("Matrix L:\n")
print(L)
print("\n\nMatrix U:\n")
print(U)


m, n = A.shape
I = np.eye(n)
# B is the inverse of A
B = np.zeros_like(A, dtype=np.float64)


# calculate columns of B
for i in range(n):
    # get the i column of I
    c = I[:,i:i+1]
    z = substitution(L, c)
    i_col_B = substitution(U, z)
    B[:,i:i+1] = i_col_B
    

A_dot_B = np.dot(B,A).astype(np.int64)    
print("\n\n  -1")
print("AA \t- Matrix multiplication of A and its inverse:")
print(A_dot_B)