import numpy as np

"""
A program that stores matrix generation functions
for various homework problems
"""

#computes a random tridiagonal matrix
def HW3P4(n):
    A = np.identity(n)
    for i in range(0,n):
        A[i,i] = np.random.random() 
        A[0,i] = np.random.random() 
        A[i,0] = np.random.random() 
    return A

#computes a central difference matrix
def HW3P5(n,val,alpha,beta):
    A = np.matrix(np.identity(n))
    h = (beta - alpha) / float((n+1))
    for i in range(0,n):
        A[i,i] = 2 + val * h * h
        if (i + 1 < n):
            A[i+1,i] = -1
            A[i,i+1] = -1

    return A

