import numpy as np
import householder

"""
This method performs the ``pure'' QR algorithm. 
Input:  B - the matrix we wish to solve
        T - the number of iterations 
"""
def pureQRalg(B,T):
    # used so we don't change our original matrix
    A = np.matrix(np.copy(B))
    # used to find the eigenvector matrix
    Qfinal = np.matrix(np.identity(A.shape[0]))
    # algorithm
    for i in range(0,T):
        Q,R = householder.householderQR(A)
        Qfinal = np.dot(Qfinal,Q)
        A = np.dot(R,Q)

    return A, Qfinal
    
