COUNT_ALL = 0
load("step1.sage")
load("Costello-Hisil.sage")
load("step2.sage")
from ctool import OpCount

# --------------- PARAMETERS ----------------------------------------------
A,p,l,k = 52, 131, 19, 1
N = 152
K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])
Q = E.random_point()

# ------------------- FIND KERNEL GENERATOR -------------------------------
# P = optimized_point_finding(A,p,k,l,N,K)
OpCount.clean()
P = point_finding(A,p,l,k)
print('P is ', P, 'order P is ', order(P))
xP = [P[0], P[2]]

# ------------------- COMPUTE COSTELLO-HISIL --------------------------------------
OpCount.clean()
kernel = kernel_points(xP, A, (l-1)/2)
ch_iso = odd_iso(kernel,Q)

# ---------------------	RUN OUR ALGO ----------------------------------------
OpCount.clean()
our_iso = evaluate_from_G(p,k,P,A,l,Q)
