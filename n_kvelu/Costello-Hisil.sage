# here we implement the Costello-Hisil algorithm 
# eprint https://eprint.iacr.org/2017/504.pdf
# for computing odd degree isogenies

def xADD(P,Q,R):
	# montgomery xADD
	xP,zP = P
	xQ,zQ = Q
	xR,zR = R
	U = (xP-zP)*(xQ+zQ)
	V = (xP+zP)*(xQ-zQ)
	res1 = zR*((U+V)**2)
	res2 = xR*((U-V)**2)
	if res2 == 0:
		res1 = 1
		res2 = 0 
	else:
		res1 = res1/res2
		res2 = 1
	return [res1, res2]


def xDBL(P,A):
	# montgomery xDBL 
	xP,zP = P
	R = (xP+zP)**2
	S = (xP-zP)**2
	t = xP*zP # T/4 from s.s. paper
	r1 = R*S
	r2 = 4*t*(S+(A+2)*t)
	if r2 == 0:
		r1 = 1
		r2 = 0 
	else:
		r1 = r1/r2
		r2 = 1
	return [r1, r2]

def criss_cross(a,b,c,d):
	# alg 1, p. 11
	# performs a small computation on the inputs
	# cost: 2M + 2a
	t1 = a*d
	t2 = b*c
	return (t1+t2, t1-t2)

def kernel_points(P, A, d): 
	# alg 2, p. 11
	# given a generator of the kernel, P, and the montgomery curve constant of the domain curve,
	# returns the first d multiples of P
	kernel = [P]
	K = parent(A)
	if d > 2:
		kernel.append(xDBL(P,A))
		for i in range(3,d):
			temp = xADD(kernel[i-2],P,kernel[i-3])
			kernel.append(temp)
	# if d==2, then kernel will have two elts: P which is added at the start, and infinity which is added at the end
	elif d == 1:
		return kernel
	kernel.append([K(1),K(0)])
	return kernel

def odd_iso(kernel,Q):
	# alg 3, p. 12
	# takes as input d kernel points, and a point to evaluate
	# returns the isogeny corresponding to input kernel, evaluated at the input point
	d = len(kernel)
	xQ,zQ = Q[0],Q[1]
	X_hat = xQ + zQ
	Z_hat = xQ - zQ
	X_prime, Z_prime = criss_cross(kernel[0][0],kernel[0][1],X_hat,Z_hat)
	for i in range(1,d):
		t0,t1 = criss_cross(kernel[i][0],kernel[i][1],X_hat,Z_hat)
		X_prime *= t0
		Z_prime *= t1
	res1 = (X_prime**2)*xQ
	res2 = (Z_prime**2)*zQ
	if res2 == 0:
		res1 = 1
		res2 = 0 
	else:
		res1 = res1/res2
		res2 = 1
	return [res1,res2]

def simultaneous_odd_iso():
	# given a set of points, and a kernel generator P, and a montgomery coefficient,
	# returns the isogeny with kernel <P> evaluated at the input points
	return 1

A,p,l,k = 5, 11, 3, 1
N = 12
K = GF(p) # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])
print('-------testing xADD---------')
P = E.random_point()
print('P is ', P)
Q = 2*P
print('Q is ', Q)
xP = [P[0], P[2]]
xQ = [Q[0], Q[2]]
R = 2*P
xR = [R[0],R[2]]
print('ours is ', xADD(xQ,xP,xP))
print('should be ', 3*P)
print('---------testing xDBL-------')
print('2P should be', 2*P)
print('ours is ', xDBL(xP,A))




