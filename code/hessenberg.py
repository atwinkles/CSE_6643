import numpy as np

"""
A program to put a matrix into Hessenberg form
Input:  B - the matrix to be transformed
"""
def hessenbergform(B):
    A = np.matrix(np.copy(B))
    N = len(A)
    for k in range(0,N-1):
        x = A[k+1:,k]
        v = np.matrix(np.sign(x[0,0])*np.linalg.norm(x)*np.eye(x.shape[0],1,0) + x)
        v = v / np.linalg.norm(v)
        A[k+1:N,k:N] = A[k+1:N,k:N] - 2 * v * (v.T * A[k+1:N,k:N])
        A[1:N,k+1:N] = A[1:N,k+1:N] - 2 * (A[1:N,k+1:N]*v) * v.T
    return A
