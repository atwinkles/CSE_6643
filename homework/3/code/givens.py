import numpy as np
import math

def givensrotation(a,b):
    r = math.hypot(a,b)
    c = a / r
    s = -b / r

    return c,s

def givensQR(A):
    m, n = A.shape
    Q = np.matrix(np.identity(m))
    R = np.matrix(np.copy(A))
    (rows, cols) = np.tril_indices(m, -1, n)
    for (row, col) in zip(rows, cols):

        # Compute Givens rotation matrix and
        # zero-out lower triangular matrix entries.
        if R[row, col] != 0:
            (c, s) = givensrotation(R[col, col], R[row, col])

            G = np.matrix(np.identity(m))
            G[[col, row], [col, row]] = c
            G[row, col] = s
            G[col, row] = -s

            R = np.dot(G, R)
            Q = np.dot(Q, G.T)

    return (Q, R)
