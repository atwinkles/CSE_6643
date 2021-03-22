import sys, time
import numpy as np
import matplotlib.pyplot as plt

import matrixgen
import iteration

mu = -5
w1 = 1
w2 = 1.5
err = 1e-10
T = 50

def graphanalysis(N,T1,T2, T3, T4 ,t1lab,t2lab,t3lab,t4lab,nval,tval,title):
    string = title + '.png'
    plt.plot(N,T1,"ro",label=t1lab)
    plt.plot(N,T2,"bo",label=t2lab)
    plt.plot(N,T3,"go",label=t3lab)
    plt.plot(N,T4,"yo",label=t4lab)
    plt.xlabel(nval)
    plt.ylabel(tval)
    plt.title(title)
    plt.legend()
    plt.savefig(string)
    plt.show() 

N = []
TJ = []
TGS = []
TSOR1 = []
TSOR2 = []

for i in range(0,5):
    n = 2**i * 100
    N.append(n)
    A = matrixgen.HW4P4(n)
    A = A + mu*np.matrix(np.identity(n))
    b = np.matrix(np.random.rand(n,1))
    x = np.matrix(np.zeros(n)).T

    start = time.time()
    sol,it = iteration.jacobi(A,x,b,err,T)
    stop = time.time() - start
    TJ.append(stop)

    start = time.time()
    sol,it = iteration.gs(A,x,b,err,T)
    stop = time.time() - start
    TGS.append(stop)

    start = time.time()
    sol,it = iteration.sor(A,x,b,w1,err,T)
    stop = time.time() - start
    TSOR1.append(stop)

    start = time.time()
    sol,it = iteration.sor(A,x,b,w2,err,T)
    stop = time.time() - start
    TSOR2.append(stop)

graphanalysis(N,TJ,TGS,TSOR1,TSOR2, "Jacobi","Gauss-Seidel","SOR, w=1","SOR, w=1.5", "Matrix Size","Average Time(s)","Iteration Comparison 1")
