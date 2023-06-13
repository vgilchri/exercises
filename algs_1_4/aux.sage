from ctool import OpCount
def xADD(P,Q,R):
	# montgomery xADD
	OpCount.op("xADD", str(k))
	xP,zP = P
	xQ,zQ = Q
	xR,zR = R
	U = (xP-zP)*(xQ+zQ)
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("mult", str(k))
	V = (xP+zP)*(xQ-zQ)
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("mult", str(k))
	res1 = zR*((U+V)**2)
	OpCount.op("square", str(k))
	OpCount.op("add", str(k))
	OpCount.op("mult", str(k))
	res2 = xR*((U-V)**2)
	OpCount.op("square", str(k))
	OpCount.op("add", str(k))
	OpCount.op("mult", str(k))
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
	OpCount.op("xDBL", str(k))
	OpCount.op("mult", str(k))
	OpCount.op("mult", str(k))
	OpCount.op("mult", str(k))
	OpCount.op("mult", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	OpCount.op("square", str(k))
	OpCount.op("square", str(k))
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
	OpCount.op("mult", str(k))
	t2 = b*c
	OpCount.op("mult", str(k))
	OpCount.op("add", str(k))
	OpCount.op("add", str(k))
	return (t1+t2, t1-t2)

#def kernel_points(P, A, d):
	# alg 2, p. 11
	# given a generator of the kernel, P, and the montgomery curve constant of the domain curve,
	# returns the first d multiples of P
#	kernel = [P]
#	if d >= 2:
		#kernel.append(xDBL(P,A))
		#for i in range(2,d):
		#	temp = xADD(kernel[i-1],P,kernel[i-2])
		#	kernel.append(temp)
	# if d==2, then kernel will have two elts: P which is added at the start, and infinity which is added at the end
	#return kernel

from ctool import OpCount

# Normal point-finding algo

def point_finding(A,p,l,k):
# checks
    if (l-1 % k == 0):
        return "k ne divise pas l-1"
    if (k > 13):
        return "k est trop grand"
    K = GF(p^k, 'x')
    print("Field: {}".format(K))
    # sample a random point
    E = EllipticCurve(K, [0, A, 0, 1, 0])
    P_l =E.random_point()
    N = E.order()
    if N % l**2 == 0:
        t = N // l**2
        OpCount.op("div", str(k))
    else:
        t = N//l
        OpCount.op("div", str(k))

    while P_l.order() != l:
        P = E.random_point()
        OpCount.op("mult", str(k))
        P_l = t*P
    return P_l

def generate_points(E, n):
    G = E.gens()[0]
    points =  []
    for i in range(1, n+1):
        P = i*G
        points.append([P[0], P[2]])
    return points
