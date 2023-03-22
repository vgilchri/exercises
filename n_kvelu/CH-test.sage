

load("step1.sage")
load("Costello-Hisil.sage")
from ctool import OpCount

# --------------- PARAMETERS ----------------------------------------------
A,p,l,k = 52, 131, 19, 1
N = 152
K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])

# ------------------- FIND KERNEL GENERATOR -------------------------------
# P = optimized_point_finding(A,p,k,l,N,K)
P = point_finding(A,p,l,k)
print('P is ', P, 'order P is ', order(P))
xP = [P[0], P[2]]

# ------------------- COMPUTE KERNEL --------------------------------------
kernel = kernel_points(xP, A, (l-1)/2)

# -------------------- FIND 2-TORSION POINTS ------------------------------
tors2 = []
for P in E:
	if order(P) == 2 and P[0] != 0:
		tors2.append(P)
print(tors2)

print('2-torsion point 1: ', odd_iso(kernel,[tors2[0][0],tors2[0][2]]))
print('2-torsion point 2: ', odd_iso(kernel,[tors2[1][0],tors2[1][2]]))

# --------------------- COMPUTE OUR CODOMAIN V.S. SAGE'S ------------------
alpha = odd_iso(kernel,[tors2[0][0],tors2[0][2]])[0]
monty_coeff = (alpha**2 + 1 )/(-1*alpha)
print('monty coeff should be ',monty_coeff)
phi = E.isogeny(P)
E1 = phi.codomain()
print('Sages codomain curve is ',E1)

