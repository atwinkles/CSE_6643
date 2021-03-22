import numpy as np
import math

"""
A program to perform QR decomposition with
Householder reflections
Input:  A - the matrix to be decomposed
"""
def householderQR(A):
    m, n = A.shape
    R = np.matrix(np.copy(A))
    Q = np.identity(m)
    for k in range(n-1):
        x = R[k:,k]
        e = np.zeros_like(x)
        e[0] = math.copysign(np.linalg.norm(x),-A[k,k])
        u = x + e
        v = u / np.linalg.norm(u)

        Qval = np.identity(m)
        Qval[k:,k:] -= 2.0 * np.outer(v,v)

        Q = np.dot(Q,Qval.T)
        R = np.dot(Qval, R)


    return Q,R
