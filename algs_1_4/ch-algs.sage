load('aux.sage')

def odd_iso(kernel,Q):
	# alg 3, p. 12
	# takes as input d kernel points, and a point to evaluate
	# returns the isogeny corresponding to input kernel, evaluated at the input point
	d = len(kernel)
	xQ,zQ = Q[0],Q[1]
	X_hat = xQ + zQ
	OpCount.op("add", str(k))
	Z_hat = xQ - zQ
	OpCount.op("add", str(k))
	X_prime, Z_prime = criss_cross(kernel[0][0],kernel[0][1],X_hat,Z_hat)
	for i in range(1,d):
		t0,t1 = criss_cross(kernel[i][0],kernel[i][1],X_hat,Z_hat)
		X_prime *= t0
		OpCount.op("mult", str(k))
		Z_prime *= t1
		OpCount.op("mult", str(k))
	res1 = (X_prime**2)*xQ
	OpCount.op("mult", str(k))
	OpCount.op("square", str(k))
	res2 = (Z_prime**2)*zQ
	OpCount.op("mult", str(k))
	OpCount.op("square", str(k))
	if res2 == 0:
		res1 = 1
		res2 = 0
	else:
		OpCount.op("inv", str(k))
		res1 = res1/res2
		res2 = 1
	return [res1,res2]

def simultaneous_odd_iso(P, A, Qs, l):
    xP = [P[0], P[2]]
    OpCount.clean()
    kernel = kernel_points(xP, A, (l-1)/2)
    kernel_hat = []
    for i in range(len(kernel)):
        OpCount.op("add", str(k))
        x_hat = kernel[i][0] + kernel[i][1]
        OpCount.op("add", str(k))
        z_hat = kernel[i][0] - kernel[i][1]
        kernel_hat.append([x_hat, z_hat])

    res = []
    for i in range(len(Qs)):
        n_p = odd_iso(kernel_hat, Qs[i])
        res.append(n_p)
    return res
