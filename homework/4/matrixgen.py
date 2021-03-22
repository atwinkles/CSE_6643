import numpy as np

"""
A program that stores matrix generation functions
for various homework problems
"""

#computes a random tridiagonal matrix
def HW4P4(m):
    A = np.matrix([[0.0 for x in range(int(m))] for y in range(int(m))])
    for i in range (m):
        A[i,i] = 10
        if i < m-1:
            A[i+1,i] = 3
            A[i,i+1] = 3
        if i < m-2:
            A[i+2,i] = 2
            A[i,i+2] = 2
        if i < m-3:
            A[i+3,i] = 1
            A[i,i+3] = 1
    return A

#computes a central difference matrix
def HW4P5(n,val,alpha,beta, ualpha, ubeta):
    A = np.matrix(np.identity(n))
    b = np.matrix(np.zeros(n)).T
    h = (beta - alpha) / float((n))
    for i in range(0,n):
        A[i,i] = 2.0 + val * h * h
        if (i + 1 < n):
            A[i+1,i] = -1.0
            A[i,i+1] = -1.0
        if (i == 0):
            b[i,0] = h*h*func(alpha + (i+1) * h) + ualpha
        elif (i == n-1):
            b[i,0] = h*h*func(alpha + (i+1) * h) + ubeta
        else:
            b[i,0] = h*h*func(alpha + (i+1) * h)

    return A, b

def func(x):
    return 3.0*x - 0.5
