import sys, time
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,'../../code/')
import matrixgen
import givens
import poweriter
import inverseiter
import pureQR

def problem4(n,t):
    N = []
    T = []
    for i in range(1,n):
        j = i * 100
        N.append(j)
        tsum = 0
        for k in range(0,t):
            print "%d, %d" % (j, k)
            A = matrixgen.HW3P4(j)
            start = time.time()
            B = givens.givensQR(A)
            tsum += (time.time() - start)
        T.append(tsum / float(t))
    graphanalysis(N,T,"Row & Column size","Average Time (s)","Given's QR Analysis")

def problem5p1(n,t):
    N = []
    T1 = []
    C1 = []
    T2 = []
    C2 = []
    for i in range(0,n):
        j = 256 * 2**i
        N.append(j)
        tsum1 = 0
        csum1 = 0
        tsum2 = 0
        csum2 = 0
        for k in range(0,t):
            A = matrixgen.HW3P5(j,0,0,1)
            start = time.time()
            maxeig = poweriter.poweriteration(A, 10)
            mineig = inverseiter.inverseiteration(A,0,10)
            tsum1 += (time.time() - start)
            csum1 += maxeig/mineig

            B = matrixgen.HW3P5(j,2,0,1)
            start = time.time()
            maxeig = poweriter.poweriteration(B, 10)
            mineig = inverseiter.inverseiteration(B,0,10)
            tsum2 += (time.time() - start)
            csum2 += maxeig/mineig
        T1.append(tsum1 / float(t))
        T2.append(tsum2 / float(t))
        C1.append(csum1 / float(t))
        C2.append(csum2 / float(t))
    graphanalysis(N,T1,T2,"lambda = 0","lambda = 2", "Matrix size","Average Time(s)", "Eigenvalue Solve Time") 
    graphanalysis(N,C1,C2,"lambda = 0","lambda = 2", "Matrix size", "Average Condition Number", "Condition Number Comparison")

def problem5p2(n,t):
    N = []
    T1 = []
    T2 = []
    for i in range(0,n):
        j = 256 * 2**i
        N.append(j)
        tsum1 = 0
        tsum2 = 0
        for k in range(0,t):
            print "%d, %d" % (j, k)
            A = matrixgen.HW3P5(j,0,0,1)
            start = time.time()
            D, Q = pureQR.pureQRalg(A, 5)
            tsum1 += (time.time() - start)

            B = matrixgen.HW3P5(j,2,0,1)
            start = time.time()
            D, Q = pureQR.pureQRalg(B, 5)
            tsum2 += (time.time() - start)
        print "done"
        T1.append(tsum1 / float(t))
        T2.append(tsum2 / float(t))
    graphanalysis(N,T1,T2,"lambda = 0","lambda = 2", "Matrix size","Average Time(s)", "QR Algorithm Time") 
        
def graphanalysis(N,T,nval,tval,title):
    plt.plot(N,T)
    plt.xlabel(nval)
    plt.ylabel(tval)
    plt.title(title)
    plt.show() 

def graphanalysis(N,T1,T2,t1lab,t2lab,nval,tval,title):
    plt.plot(N,T1,"o",label=t1lab)
    plt.plot(N,T2,"o",label=t2lab)
    plt.xlabel(nval)
    plt.ylabel(tval)
    plt.title(title)
    plt.legend()
    plt.show() 
