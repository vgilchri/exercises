import sys
import random
from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_odd

from ctool import OpCount

def create_elliptic_curve(k):
    while True:
      a1 = 0
      a2 = k.random_element()
      a3 = 0
      a4 = k.random_element()
      a6 = k.random_element()
      try:
        E = EllipticCurve([a1,a2,a3,a4,a6])
        if E.cardinality().is_power_of(2):
            continue
        else:
            return E, [a1,a2,a3,a4,a6]
      except:
        pass

def find_point_of_order(E):
     while True:
       P = E.random_point()
       l = P.order()
       # paper states velu formulas for odd prime l
       # but this script allows any odd order l (including 1)
       if l%2:
           return P

def compute_rationals(E, P, kX, k, a_inv):
    a1,a2,a3,a4,a6 = a_inv
    l = P.order()
    x = {s:(s*P)[0] for s in range(1,l)}

    S = range(1,(l+1)/2)
    if l > 1:
      assert max(S) == (l-1)/2

    YY = X^3+a2*X^2+a4*X+a6
    YYprime = YY.diff()

    Psi = prod(X-x[s] for s in S)
    Psi = kX(Psi) # in case l=1
    Psiprime = Psi.diff()
    Psiprimeprime = Psiprime.diff()

    xsum = sum(x[s] for s in range(1,l))
    Phi = 4*YY*(Psiprime^2-Psiprimeprime*Psi)-2*YYprime*Psiprime*Psi+(l*X-xsum)*Psi^2
    assert Phi == 4*(X^3+a2*X^2+a4*X+a6)*(Psiprime^2-Psiprimeprime*Psi)-2*(3*X^2+2*a2*X+a4)*Psiprime*Psi+(l*X-xsum)*Psi^2
    Phiprime = Phi.diff()

    Omega = Phiprime*Psi-2*Phi*Psiprime

    print('try','velu',l,E,P,Psi,Phi,Omega)
    return Psi, Phi, Omega

def build_new_curve(kX, Psi, Phi, Omega, E, P, a_inv):
    a1,a2,a3,a4,a6 = a_inv
    sageisog = E.isogeny(P)
    l = P.order()
    print("Psi: {}".format(Psi))
    print("\n Kpoly: {}".format(sageisog.kernel_polynomial()(X)))
    assert sageisog.kernel_polynomial()(X) == Psi
    S = range(1,(l+1)/2)
    x = {s:(s*P)[0] for s in range(1,l)}
    xsum = sum(x[s] for s in range(1,l))

    s1 = -Psi[len(S)-1]
    assert s1 == sum(x[s] for s in S)
    assert 2*s1 == xsum
    s2 = Psi[len(S)-2]
    s3 = -Psi[len(S)-3]

    b2 = 4*a2
    b4 = 2*a4
    b6 = 4*a6
    assert E.b_invariants()[:3] == (b2,b4,b6)

    t = 6*(s1^2-2*s2)+b2*s1+b4*len(S)
    w = 10*(s1^3-3*s1*s2+3*s3)+2*b2*(s1^2-2*s2)+3*b4*s1+b6*len(S)

    assert (t,w) == compute_vw_kohel_odd(b2,b4,b6,s1,s2,s3,len(S))

    newa1 = 0
    newa2 = a2
    newa3 = 0
    newa4 = a4-5*t
    newa6 = a6-b2*t-7*w

    assert sageisog.codomain().a_invariants() == (newa1,newa2,newa3,newa4,newa6)
    new_as = [newa1,newa2,newa3,newa4,newa6]
    E_phi = EllipticCurve([newa1,newa2,newa3,newa4,newa6])
    return E_phi, new_as

def sanity_check(E, a_inv, E_phi, New_as,  Phi, Omega, Psi, k, kX):
    a1,a2,a3,a4,a6 = a_inv
    newa1,newa2,newa3,newa4,newa6 = New_as
    newX = Phi/Psi^2
    newYY = YY*Omega^2/Psi^6
    assert newYY == newX^3+newa2*newX^2+newa4*newX+newa6

     # cyclotomic statement

    for zeta in k:
      if zeta == 0: continue
      n = zeta.multiplicative_order()
      h = prod(X-zeta^i for i in range(n) if gcd(i,n) == 1)
      h = kX(h)
      assert h == cyclotomic_polynomial(n)(X)
