import time
import numpy as np
from math import hypot

def toepmatrixgen(m):
    A = np.array([[0.0 for x in range(int(m))] for y in range(int(m))])
    for i in range (m):
        A[i][i] = 10
        if i < m-1:
            A[i+1][i] = 3
            A[i][i+1] = 3
        if i < m-2:
            A[i+2][i] = 2
            A[i][i+2] = 2
        if i < m-3:
            A[i+3][i] = 1
            A[i][i+3] = 1
    return A

def svdmatrixgen(m):
    H = np.random.randn(m,m)
    Q = modifiedGS(H)

    A = np.array([[0.0 for x in range(int(m))] for y in range(int(m))])
    for i in range(m):
        A[i][i] = 3.0**-(i+1)

    A = np.matmul(Q,(np.matmul(A,np.matrix.transpose(Q))))

    return A

def twonorm(A):
    val = 0
    for i in A:
        val += i*i
    return np.sqrt(val) 

def compare(A,B,q):
    val = 0

    if q == 1:
        for i in range(len(A)):
            for j in range(len(A)):
                val += abs(A[i][j] - B[i][j])

    if q == 2:
        for i in range(len(A)):
            for j in range(len(A)):
                val += (A[i][j] - B[i][j])**2  


    return val

def classicGS(A):
    m = len(A)
    Q = np.array([[0.0 for x in range(int(m))] for y in range(int(m))])
    R = np.array([[0.0 for x in range(int(m))] for y in range(int(m))])
    for j in range(m):
        v = A[j]
        for i in range(j-1):
            R[i][j] = Q[i].T.dot(A[j])
            v = v - R[i][j]*Q[i]
        R[j][j] = twonorm(v)
        Q[j] = v/R[j][j]
    return Q

def modifiedGS(A):
    m = len(A)
    Q = np.array([[0.0 for x in range(int(m))] for y in range(int(m))])
    R = np.array([[0.0 for x in range(int(m))] for y in range(int(m))])
    V = np.array([[0.0 for x in range(int(m))] for y in range(int(m))])
    for i in range(m):
        V[i] = A[i]
    for i in range(m):
        R[i][i] = twonorm(V[i])
        Q[i] = V[i] / R[i][i]
        for j in range(i+1,m):
            R[i][j] = Q[i].T.dot(V[j])
            V[j] = V[j] - R[i][j]*Q[i]
    return Q

def householderQR(A):
    m = len(A)
    Q = np.array([[0.0 for x in range(int(m))] for y in range(int(m))])
    for k in range(m):
        for l in range(m):
            Q[k][l] = A[k][l]
    for k in range(m):
        x = Q[k:m,k]
        v = np.sign(x[0]) * twonorm(x)*np.eye(1,m-k) + x
        v = v / twonorm(v)
        Q[k:m,k:m] = Q[k:m,k:m] - 2*v*(v.T*Q[k:m,k:m])
    return Q
        
def givensQR(A):
    m = len(A)
    Q = np.identity(m)
    R = np.copy(A)
    for i in range(m):
        for j in range(m):
            if R[i][j] != 0:
                r = hypot(R[j][j],R[i][j])
                c = R[j][j]/r
                s = - R[i][j] / r
                G = np.identity(m)
                G[[j,i],[j,i]] = c
                G[i][j] = s
                G[j][i] = -s

                R = np.dot(G,R)
                Q = np.dot(Q,G.T)
    return Q
