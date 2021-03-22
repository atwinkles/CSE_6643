import numpy as np

def centraldiff(n,val,alpha,beta,ualpha,ubeta):
    A = np.matrix(np.identity(n))
    B = np.matrix(np.zeros((n,1)))
    h = (beta - alpha) / float((n+1))
    for i in range(0,n):
        A[i,i] = 2 + val * h * h
        if (i + 1 < n):
            A[i+1,i] = -1
            A[i,i+1] = -1
        if (i == 1):
            B[i] = h * h * func(alpha + i*h) + ualpha
        elif (i == n):
            B[i] = h * h * func(alpha + i*h) + ubeta
        else:
            B[i] = h * h * func(alpha + i*h)

    return A, B

def func(x):
    return 3.0 * x - 1.0/2.0
