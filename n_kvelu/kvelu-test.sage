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

# ---------------------	RUN OUR ALGO ----------------------------------------



for line in open(filename, 'r'):
    line = line.split(',')
    k = int(line[0].replace("[", ""))
    A = int(line[1])
    p = int(line[2])
    l = int(line[3].replace("]", ""))

    if k%2 != 0 and k < 12:
        print("A, p, l, k, : {} {} {} {} ".format(A, p, l, k))
        K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
        A = GF(p)(A)
        E = EllipticCurve(K, [0, A, 0, 1, 0])
        Q = E.random_point()

        P = point_finding(A,p,l,k)
        OpCount.clean()
        our_iso = evaluate_from_G(p,k,P,A,l,Q)
        print("A, p, l, k, Q, iso: {} {} {} {} {} {}".format(A, p, l, k, Q, our_iso))
