import numpy as np

"""
Uses the power method to find the largest
eigenvalue of a matrix.
Inputs: A - the matrix to be solved
        T - the number of iterations
"""
def poweriteration(A, T):
    N = A.shape[0]
    v = np.matrix(np.eye(N,1,0))
    for i in range(1,T):
        w = A * v
        v = w / np.linalg.norm(w)
        eig = float(v.T * A * v)

    return eig
