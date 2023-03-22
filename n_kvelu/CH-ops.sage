load("step1.sage")
load("Costello-Hisil.sage")
from ctool import OpCount

A,p,l,k = 52, 131, 19, 1
N = 152
K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])

# ---------- FIND KERNEL GENERATOR ------------------
# P = optimized_point_finding(A,p,k,l,N,K)
P = point_finding(A,p,l,k)
print('P is ', P, ', order P is ', order(P))
xP = [P[0], P[2]]

# ----------- COMPUTE KERNEL -----------------------
OpCount.clean() 
kernel = kernel_points(xP, A, (l-1)/2)
print('kernel is ', kernel)
OpCount.print_results()
OpCount.clean() 

# ------------- EVALUATE AT RANDOM POINT -----------

Q = E.random_point()
xQ = [Q[0],Q[2]]
print('Q is ', xQ)
phi_Q = odd_iso(kernel,xQ)
print('phi(Q) is ',phi_Q)
OpCount.print_results()
OpCount.clean() 