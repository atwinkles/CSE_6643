import numpy as np

"""
Uses the inverse method to find the largest
eigenvalue of a matrix.
Inputs: A - the matrix to be solved
        mu - the eigenvalue guess
        T - the number of iterations
"""
def inverseiteration(A,mu,T):
    N = A.shape[0]
    v = np.matrix(np.eye(N,1,0))
    B = np.linalg.inv(A - mu * np.identity(N))
    for i in range(1,T):
        w = B * v
        v = w / np.linalg.norm(w)
        eig = float(v.T * A * v)

    return eig
