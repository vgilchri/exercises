COUNT_ALL = 1
load("step1.sage")
load("step2.sage")
from ctool import OpCount

# --------------- PARAMETERS ----------------------------------------------
A,p,l,k = 52, 131, 19, 1
N = 152


K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])
filename = "params.txt"
for line in open(filename, 'r'):
    line = line.split(',')
    A = int(line[0].replace("(", ""))
    p = int(line[1])
    l = int(line[2])
    k = int(line[3])
    if k%2 != 0 and k < 12:
        j_inv = int(line[4])
        ratio = int(line[5].replace(")", ""))
        K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
        A = GF(p)(A)
        E = EllipticCurve(K, [0, A, 0, 1, 0])
        Q = point_finding(A,p,l,k)
        print("A, p, l, k, j_inv, ratio, Q: {} {} {} {} {} {} {}".format(A, p, l, k, j_inv, ratio, Q))
