
load("step1.sage")
load("Costello-Hisil.sage")

A,p,l,k = 5, 11, 2, 1
N = 12
K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])

# P = optimized_point_finding(A,p,k,l,N,K)
P = point_finding(A,p,l,k)
print('P is ', P, 'order P is ', order(P))
xP = [P[0], P[2]]
kernel = kernel_points(xP, A, l)
print('our kernel is ', kernel)
their_kernel = []
for i in range(1,l+1):
	their_kernel.append(i*P)
print('their kernel is ', their_kernel)

Q = E.random_point()
xQ = [Q[0],Q[2]]
#print(parent(xQ), 'and ', parent(xP))
print('Q is ', Q)
new_Q = odd_iso(kernel, xQ)
print('our phi(Q) is ', new_Q)

phi = E.isogeny(P)
their_Q = phi(Q)
print('their phi(Q) is ', their_Q)
E1 = phi.codomain()

for R in E1:
	print(R, ' has order ', order(R))



