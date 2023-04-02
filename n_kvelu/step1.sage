#*****************************************************************************

#       Algorithm for step 1: find a point G of order \ell

#*****************************************************************************

from ctool import OpCount

# for comparison we start with the normal point-finding algo

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
    t = E.order()//l
    while P_l.order() != l:
        P = E.random_point()
        P_l = t*P
    return P_l

# computes a point evaluated at the frobenius map

def frob_power(K,E,P,power):
    frob = K.frobenius_endomorphism(power)
    OpCount.op("frob", str(K.degree()))
    OpCount.op("frob", str(K.degree()))
    R = E(frob(P[0]), frob(P[1]))
    return R

# returns the image of a point P in E by the endo delta

# delta is defined as the division of (X^k - 1) and the kth cyclotomic polynomial

def delta(P,E,k,p):
    R = "error, wrong k"
    K = GF(p**k, 'x')
    if k == 1:
        return P
    if (k in {2,3,5,7,11}):
        R = frob_power(K,E,P,1) - P
    elif k==4 :
        R = frob_power(K,E,P,2) - P
    elif k==6 :
        S = frob_power(K,E,P,3) - P
        R = frob_power(K,E,S,1) + S
    elif k==8 :
        R = frob_power(K,E,P,4) - P
    elif k==9 :
        R = frob_power(K,E,P,3) - P
    elif k==10 :
        S = frob_power(K,E,P,5) - P
        R = frob_power(K,E,S,1) + S
    elif k==12 :
        S = frob_power(K,E,P,6) - P
        R = frob_power(K,E,S,2) + S
    return R

def get_h_k(A,p,k) :

    K = GF(p**k)

    E = EllipticCurve(K, [0, A, 0, 1, 0])

    Nk = E.order()
    N = Nk

    if (k in {2,3,5,7,11}) :

        K1 = GF(p)

        E1 = EllipticCurve(K1, [0, A, 0, 1, 0])

        N1 = E1.order()

        N= Nk//N1

    elif k==4 :

        K2 = GF(p**2)

        E2 = EllipticCurve(K2, [0, A, 0, 1, 0])

        N2 = E2.order()

        N= Nk//N2

    elif k==6 :

        K1 = GF(p)

        E1 = EllipticCurve(K1, [0, A, 0, 1, 0])

        N1 = E1.order()

        K2 = GF(p**2)

        E2 = EllipticCurve(K2, [0, A, 0, 1, 0])

        N2 = E2.order()

        K3 = GF(p**3)

        E3 = EllipticCurve(K3, [0, A, 0, 1, 0])

        N3 = E3.order()

        N= (Nk*N1)//(N2*N3)

    elif k==8 :

        K4 = GF(p**4)

        E4 = EllipticCurve(K4, [0, A, 0, 1, 0])

        N4 = E4.order()

        N= Nk//N4

    elif k==9 :

        K3 = GF(p**3)

        E3 = EllipticCurve(K3, [0, A, 0, 1, 0])

        N3 = E3.order()

        N= Nk//N3

    elif k==10 :

        K1 = GF(p)

        E1 = EllipticCurve(K1, [0, A, 0, 1, 0])

        N1 = E1.order()

        K2 = GF(p**2)

        E2 = EllipticCurve(K2, [0, A, 0, 1, 0])

        N2 = E2.order()

        K5 = GF(p**5)

        E5 = EllipticCurve(K5, [0, A, 0, 1, 0])

        N5 = E5.order()

        N= (Nk*N1)//(N2*N5)

    elif k==12 :

        K2 = GF(p**2)

        E2 = EllipticCurve(K2, [0, A, 0, 1, 0])

        N2 = E2.order()

        K4 = GF(p**4)

        E4 = EllipticCurve(K4, [0, A, 0, 1, 0])

        N4 = E4.order()

        K6 = GF(p**6)

        E6 = EllipticCurve(K6, [0, A, 0, 1, 0])

        N6 = E6.order()

        N= (Nk*N2)//(N4*N6)

    return N

# input: (montgomery coefficient A for curve E) and (order l) and (finite field p and k) and

# (number of points Nratio=N/N1) and (m=Nratio/l), output: (G a point of order l on E in K)

# assume inputs are well defined for the context

def optimized_point_finding(A,p,k,l,N, K):

    # checks

    if l-1 % k == 0:
        return "error: k ne divise pas l-1"

    if k > 13:
        return "error: k est trop grand"
    # sample a random point

    E = EllipticCurve(K, [0, A, 0, 1, 0])

    P_l = E(0)


    t = N//l


    while P_l.order() != l:
        P = E.random_point()
        P = delta(P,E,k,p)
        P_l = t*P
        OpCount.op("mult", str(k))


    return P_l
