#*****************************************************************************

#       Algorithm for step 1: find a point G of order \ell

#*****************************************************************************

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
    OpCount.op("rand", str(k))
    N = E.order()
    if N % l**2 == 0:
        t = N // l**2
        OpCount.op("div", str(k))
    else:
        t = N//l
        OpCount.op("div", str(k))

    while P_l.order() != l:
        P = E.random_point()
        OpCount.op("rand", str(k))
        OpCount.op("mult", str(k))
        P_l = t*P
    return P_l


