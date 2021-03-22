import numpy as np

"""
This file holds iteration methods for solving Ax = b
"""

"""
The Jacobi method solves Ax = b
Inputs:
    A   -   The matrix
    x   -   The initial guess vector
    b   -   The vector
    err -   The error desired for the solution
    T   -   The number of iterations to be done
"""
def jacobi(A,x,b,err,T):
    n = A.shape[0]
    count = 0
    while (count < T+1):
        if count % 10 == 0:
            print "%d \t %f \t %f" % (count, x[0],x[1])
        u = np.matrix(np.zeros(n)).T
        for i in range(n):
            s = 0.0
            for j in range(n):
                if j != i:
                    s += float(A[i,j]*x[j])
            u[i] = b[i] - s
            u[i] /= float(A[i,i])

        x = u
        if (np.linalg.norm(b - A*x) < err):
            return x
        count += 1

    return x,count
            
def gs(A,x,b,err,T):
    n = A.shape[0]
    count = 0
    x = x.astype(float)
    while (count < T+1):
        if count % 10 == 0:
            print "%d \t %f \t %f" % (count, x[0],x[1])
        for i in range(n):
            s = 0.0
            for j in range(n):
                if j != i:
                    s += float(A[i,j])*float(x[j,0])
            k = (b[i,0] - s) / float(A[i,i])
            x[i,0] = k

        if (np.linalg.norm(b - A*x) < err):
            return x
        count += 1

    return x

def sor(A,x,b,w,err,T):
    n = A.shape[0]
    count = 0
    x = x.astype(float)
    while (count < T+1):
        if count % 10 == 0:
            print "%d \t %f \t %f" % (count, x[0],x[1])
        for i in range(n):
            s = 0.0
            for j in range(n):
                if j != i:
                    s += float(A[i,j]*x[j])
            x[i] = (1-w)*x[i] + w*float(b[i] - s)/float(A[i,i])

        if (np.linalg.norm(b - A*x) < err):
            return x
        count += 1

    return x
     
def steepestdescent(A,x,b,err,T):
    count = 0
    x = x.astype(float)
    while (count < T+1):
        if count % 10 == 0:
            print "%d \t %f \t %f" % (count, x[0],x[1])
        r = b - A*x
        alpha = float(r.T*r / (r.T*A*r))
        x = x + alpha * r
        if (np.linalg.norm(b - A*x) < err):
            return x
        count += 1
    return x

def gradiantdescent(A,x,b,err,T):
    count = 0
    x = x.astype(float)
    r = b - A*x
    p = r
    while (count < T+1):
        if count % 10 == 0:
            print "%d \t %f \t %f" % (count, x[0],x[1])
        alpha = float(r.T*r / (p.T*A*p))
        x = x + alpha * p
        t = 1/float(r.T*r)
        r = r - alpha * A * p
        p = r + float(r.T*r) * t * p
        if (np.linalg.norm(b - A*x) < err):
            return x
        count += 1
    return x
