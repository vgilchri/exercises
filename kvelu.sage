import sys

import random

from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_vw_kohel_odd



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

    return Psi, Phi, Omega

def build_new_curve(kX, Psi, Phi, Omega, E, P, a_inv):

    a1,a2,a3,a4,a6 = a_inv

    sageisog = E.isogeny(P)

    l = P.order()

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





#*****************************************************************************

#       Algorithm for step 1: find a point G of order \ell

#*****************************************************************************



# for comparison we start with the normal point-finding algo

def point_finding(A,p,l,k):

    # checks

    if (l-1 % k == 0):

        return "k ne divise pas l-1"

    if (k > 13):

        return "k est trop grand"



    K.<x> = GF(p, k)

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

    R = E(frob(P[0]), frob(P[1]))

    return R

# returns the image of a point P in E by the endo delta

# delta is defined as the division of (X^k - 1) and the kth cyclotomic polynomial

def delta(P,E,k,p):

    R = "error, wrong k"

    K = GF(p**k)

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

def optimized_point_finding(A,p,k,l,N):

    # checks

    if l-1 % k == 0:

        return "error: k ne divise pas l-1"

    if k > 13:

        return "error: k est trop grand"



    K = GF(p**k) # finite field of size p^k



    # sample a random point

    E = EllipticCurve(K, [0, A, 0, 1, 0])

    P_l =E(0)

    t = N//l

    while P_l.order() != l:
    	print("ok")
    	P = E.random_point()
    	P = delta(P,E,k,p)
    	P_l = t*P

    return P_l

# test



#*****************************************************************************

#       Algorithm for Step 2 : Computing an \ell-Isogeny

#*****************************************************************************



#####################################################

########## Valuation of an integer n wrt 2 ##########

#####################################################



def valuation_2(n) :

  res=0

  while n%2==0 :

    res+=1

    n=n/2

  return (res)

#####################################

########## Montgomery xADD ##########

#####################################



def xADD(xP,xQ,xR):

    return (xP*xQ-1)^2/(xP-xQ)^2/xR

#####################################

########## Montgomery xDBL ##########

#####################################



def xDBL(x,a):

    return (x^2-1)^2/(4*x*(x^2+a*x+1))


#####################################

############## Lemma 3 ##############

#####################################





def Condition_lemma3(k_prime,k_2,k_2_prime,l) :

  b= False

  if k_2==l-1 :

    b=True

  if k_2==(l-1)/2 :

    if l%4==3 :

      b=True

    elif valuation_2(k_prime)==valuation_2(k_2_prime)+1 :

      b==True



  return(b)

#####################################

########### Lemma 3 primes ##########

#####################################



def list_primes_lemma(n):

  list_of_fields = []

  for l in primes(100000):

      K = GF(l)

      element = K.2

      if element.multiplicative_order() == l-1:

        list_of_fields.append(l)



  return(list_of_fields)

#       Version where we get the full kernel polynomial and then evaluate

#*****************************************************************************



def get_S0_true(G,A,l,k_prime):

  S_0=[G[0]]

  x=G[0]

  for i in range (1,(l-1)/(2*k_prime)) :

    x=xDBL(x,A)

    S_0.append(x)



  return S_0



def get_S0_T_false(G,A,l,k_prime,K,E):

  S_0=[G[0]]

  T=[]

  x0=G

  frob=K.frobenius_endomorphism()

  x1=G

  for i in range (k) :

    T.append(x1[0])

    x1=E(frob(x1[0]), frob(x1[1]))



  while len(S_O) < (l-1)/2 :

    while x0[0] in T :

      x0=x0+G

    x1=x0

    S_0.append(x0[0])

    for i in range (k) :

        T.append(x1[0])

        x1=E(frob(x1[0]), frob(x1[1]))



  return T, S_0

#####################################

######### Get kernel points #########

#####################################



def get_kernel_polynomial_points (p,K,G,A,l, k) :

  E = EllipticCurve(K,[0,A,0,1,0])

  k_prime=p

  q = pow(p,k)

  if k%2==0 :

    k_prime=k/2

  K2=GF(l)

  element = K2(2)

  k_2=element.multiplicative_order()

  k_2_prime=k_2

  if k_2%2==0 :

    k_2_prime=k_2/2

  print(Condition_lemma3(k_prime,k_2,k_2_prime,l))

  if Condition_lemma3(k_prime,k_2,k_2_prime,l) :

    S_0=get_S0_true(G,A,l,k_prime)

    T=copy(S_0)

    for i in range (1,k) :
    	for j in range (len(S_0)) :
    		T.append(pow(S_0[j],q))
  else :
    T, S_0=get_S0_T_false(G,A,l,k_prime,K,E)

  return(T)



#def Evaluate_from_T(T,K,P) :
#  #P=K(P)
#  res=P[0]-T[0]
#  for i in  range (1,len(T)) :
#  	print(res)
# 	res= res*(P[0]-T[i])
#
#  return(res)




#       Version where we get evaluate without computing the whole kernel polynomial

#****************************************************************************************



def get_S0_false(G,A,l,k_prime,K,E):

  S_0=[G]

  lam= get_frob_value()

  r=1

  T_L=[[0,0,0] for i in range (l-1)] #[b,m,n] where b is 1 if it has been encountered, m i is the Galois orbit generator, n is its frobenius powering

  for i in range (k) :

    power=pow(lam,i)%l

    T_L[power-1]=[1,0,i]



  while len(S_0)< (l-1)/k :

    while T_L[r-1][0]==1 :

      r+=1



    previous=S_0[T_L[r-2][1]]*(pow(lam,T_L[r-2][2])%l)

    previous_2=S_0[T_L[r-3][1]]*(pow(lam,T_L[r-3][2])%l)

    new= xADD(previous,G,previous_2)

    S_0.append(new)

    for i in range (k) :

      power=pow(lam,i)%l

      T_L[(r*power)%l-1]=[1,r,i]


    return(S_0)



def evaluate_from_G(p,k,G,A,l,P):

    K = GF((p,k),'x')

    E = EllipticCurve(K,[0,A,0,1,0])

    k_prime=k

    if k%2==0 :

    	k_prime=k/2

    K2=GF(l)

    element = K2(2)

    k_2=element.multiplicative_order()

    k_2_prime=k_2

    if k_2%2==0 :

    	k_2_prime=k_2/2

    if Condition_lemma3(k_prime,k_2,k_2_prime,l) :

    	S_0=get_S0_true(G,A,l,k_prime)

    else :

    	S_0=get_S0_false()

    res=P[0]-S_0[0]

    for i in range (1,len(S_0)) :

    	res=res*(P[0]-S_0[i])

    power=p

    c=p

    for i in range (2,k) :

    	c=c*p

    	power=power + c

    res=pow(res,power)

    return(res)


#       Get the isogeny

#****************************************************************************************

def compute_rationals_bis(E, T, l, kX, k, A):

    a1,a2,a3,a4,a6 = 0,A,0,1,0

    S = range(1,(l+1)/2)

    if l > 1:

      assert max(S) == (l-1)/2



    YY = X^3+a2*X^2+a4*X+a6

    YYprime = YY.diff()



    Psi = prod(X-t for t in T)

    Psi = kX(Psi) # in case l=1

    Psiprime = Psi.diff()

    Psiprimeprime = Psiprime.diff()



    xsum = sum(t for t in T)

    Phi = 4*YY*(Psiprime^2-Psiprimeprime*Psi)-2*YYprime*Psiprime*Psi+(l*X-xsum)*Psi^2

    assert Phi == 4*(X^3+a2*X^2+a4*X+a6)*(Psiprime^2-Psiprimeprime*Psi)-2*(3*X^2+2*a2*X+a4)*Psiprime*Psi+(l*X-xsum)*Psi^2

    Phiprime = Phi.diff()



    Omega = Phiprime*Psi-2*Phi*Psiprime

    return Psi, Phi, Omega


# voir dans le cas k pair si ya pas des trucs à gérer


print("test for elliptic curves constructions")
K = GF((17,4),'x')

q = 17
k.<fieldgen> = GF(q)

kX.<X> = k[]

E, a_inv = create_elliptic_curve(k)

P = find_point_of_order(E)

print("EllipticCurve: {}\nkX: {}\nP: {}".format(E, kX, P))

Psi, Phi, Omega = compute_rationals(E, P, kX, k, a_inv)

print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))

E_phi, new_as = build_new_curve(kX, Psi, Phi, Omega, E, P, a_inv)

print("New EllipticCurve: {}".format(E_phi))


#------ SAGE CODE #-------

sageisog = E.isogeny(P)

print("sageiso: {}".format(sageisog))

A= 8
l= 5
k= 4
p=17

# test

print("test for point finding")

P=point_finding(A, p, l, k)

K = GF((p,k),'x')

E = EllipticCurve(K,[0,A,0,1,0])


print(P)

print(P.order())

print("test for optimized point finding")

N=get_h_k(A, p, k)
print(N)

E = EllipticCurve(K,[0,A,0,1,0])

P1=optimized_point_finding(A, p, k,l,N)

#print(E.cardinality()/5)

print(P1)

print(P1.order())

print("test for get_kernel_polynomial_points")

T= get_kernel_polynomial_points(p,K,P,A,l, k)
print(T)
print(len(T))
print("T[0]=",T[0])

Prand=E.random_point()
print("Prand[0]=",Prand[0])
#poly=Evaluate_from_T(T,K,Prand)
#print("?")
#print(poly)



print("test for elliptic curves constructions bis")



print("EllipticCurve: {}\nkX: {}\nP: {}".format(E, kX, P))

Psi, Phi, Omega = compute_rationals_bis(E, T, l, kX, k, A)

print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))

a_inv= 0,A,0,1,0
Psi, Phi, Omega = compute_rationals(E, P, kX, k, a_inv)

print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))


E_phi, new_as = build_new_curve(kX, Psi, Phi, Omega, E, P, a_inv)

print("New EllipticCurve: {}".format(E_phi))

