import QR
import time

for i in range(10,500,10):
    A = QR.svdmatrixgen(i)
    start = time.time()
    Q1 = QR.householderQR(A)  
    Q2 = QR.modifiedGS(A)
    diff1 = QR.compare(Q1,Q2,1)
    diff2 = QR.compare(Q1,Q2,2)
    stop = (time.time() - start)

    print "For size %d * %d, the computation took %f time, with differences norm 1: %d and norm 2: %d" % (i, i, stop, diff1, diff2)
