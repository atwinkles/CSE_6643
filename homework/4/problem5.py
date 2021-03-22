import sys, time
import numpy as np
import matplotlib.pyplot as plt
from math import exp, sqrt

import matrixgen
import iteration

err = 1e-5
T = 2000

def graphanalysis(N,T1,T2 ,t1lab,t2lab,nval,tval,title):
    string = title + '.png'
    plt.plot(N,T1,"ro",label=t1lab)
    plt.plot(N,T2,"bo",label=t2lab)
    plt.xlabel(nval)
    plt.ylabel(tval)
    plt.title(title)
    plt.legend(loc="left")
    plt.savefig(string)
    plt.show()

def graphanalysis2(N,T1,T2 ,t1lab,t2lab,nval,tval,title):
    string = title + '.png'
    plt.plot(N,T1,"r",label=t1lab)
    plt.plot(N,T2,"b",label=t2lab)
    plt.xlabel(nval)
    plt.ylabel(tval)
    plt.title(title)
    plt.legend(loc="left")
    plt.savefig(string)
    plt.show()

N = []
TSD0 = []
TCG0 = []

TSD2 = []
TCG2 = []

for i in range(5,6):
    n = 2**i * 100
    N.append(n)
    A, b = matrixgen.HW4P5(n, 0, 0, 1, 0, -2)
    x = np.matrix(np.zeros(n)).T

    start = time.time()
    sol1, it = iteration.steepestdescent(A,x,b,err,T)
    stop = time.time() - start
    TSD0.append(stop)

    start = time.time()
    sol2, it = iteration.gradiantdescent(A,x,b,err,T)
    stop = time.time() - start
    TCG0.append(stop)

    A, b = matrixgen.HW4P5(n, 2, 0, 1, 0, -2)

    start = time.time()
    sol3, it = iteration.steepestdescent(A,x,b,err,T)
    stop = time.time() - start
    TSD2.append(stop)

    start = time.time()
    sol4, it = iteration.gradiantdescent(A,x,b,err,T)
    stop = time.time() - start
    TCG2.append(stop)

def val0(i):
    return 1.0/4.0*i*(-2*i*i + i - 7)

def val2(i):
    return (exp(-sqrt(2)*i)*(exp(sqrt(2)*i)*(1-6*i)-exp(2*sqrt(2)*i)-13*exp(sqrt(2)*(2*i+1)) + exp(sqrt(2)*(i+2))*(6*i-1)+exp(2*sqrt(2))+13*exp(sqrt(2))))/(4*(exp(2*sqrt(2))-1))

T = []

for i in range(0,len(sol1)):
    T.append(0 + i * (1/float(len(sol1)+1)))

V0 = []
V2 = []

print sol1

for i in T:
    V0.append(val0(i))
    V2.append(val2(i))

print sol1

graphanalysis(N,TSD0, TCG0, "Steepest Descent", "Conjugate Gradient", "Matrix Size", "Time (s)", "Time Comparison for lambda = 0")
graphanalysis(N,TSD2, TCG2, "Steepest Descent", "Conjugate Gradient", "Matrix Size", "Time (s)", "Time Comparison for lambda = 2")

graphanalysis2(T,V0, sol1, "Exact Solution", "Approximate Solution", "x", "u", "Steepest Descent Comparison for lambda = 0")
graphanalysis2(T,V2, sol3, "Exact Solution", "Approximate Solution", "x", "u", "Steepest Descent Comparison for lambda = 2")

graphanalysis2(T,V0, sol2, "Exact Solution", "Approximate Solution", "x", "u", "Conjugate Gradient Comparison for lambda = 0")
graphanalysis2(T,V2, sol4, "Exact Solution", "Approximate Solution", "x", "u", "Conjugate Gradient Comparison for lambda = 2")
