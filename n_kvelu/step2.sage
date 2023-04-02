#*****************************************************************************

#       Algorithm for Step 2 : Computing an \ell-Isogeny

#*****************************************************************************

from ctool import OpCount
# we count xADD's, xDBL's, add's, mult's, div's, square's

#####################################################

########## Valuation of an integer n wrt 2 ##########

#####################################################



def valuation_2(n) :
  res=0
  while (n%2)==0 :
    res+=1
    n=n/2
  return (res)

#####################################

########## Montgomery xADD ##########

#####################################



#def xADD(xP,xQ,xR):
#    return (xP*xQ-1)^2/(xP-xQ)^2/xR

def xADD(P,Q,R):
	# montgomery xADD
    OpCount.op("xADD", str(k))
    xP,zP = P
    xQ,zQ = Q
    xR,zR = R
    U = (xP-zP)*(xQ+zQ)
    if COUNT_ALL:
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        OpCount.op("mult", str(k))
        V = (xP+zP)*(xQ-zQ)
    if COUNT_ALL:
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        OpCount.op("mult", str(k))
    res1 = zR*((U+V)**2)
    if COUNT_ALL:
        OpCount.op("square", str(k))
        OpCount.op("add", str(k))
        OpCount.op("mult", str(k))
        res2 = xR*((U-V)**2)
    if COUNT_ALL:
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


#####################################

########## Montgomery xDBL ##########

#####################################



def xDBL(x,a):
    if COUNT_ALL:
        OpCount.op("mult", str(k))
        OpCount.op("mult", str(k))
        OpCount.op("mult", str(k))
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        OpCount.op("add", str(k))
        OpCount.op("div", str(k))
        OpCount.op("square", str(k))
        OpCount.op("square", str(k))
        OpCount.op("square", str(k))
    OpCount.op("xDBL", str(k))
    return (x^2-1)^2/(4*x*(x^2+a*x+1))



#def xDBL(P,A):
#	# montgomery xDBL
#	xP,zP = P
#	R = (xP+zP)**2
#	S = (xP-zP)**2
#	t = xP*zP # T/4 from s.s. paper
#	r1 = R*S
#	r2 = 4*t*(S+(A+2)*t)
#	OpCount.op("mult", str(k))
#	OpCount.op("mult", str(k))
#	OpCount.op("mult", str(k))
#	OpCount.op("mult", str(k))
#	OpCount.op("add", str(k))
#	OpCount.op("add", str(k))
#	OpCount.op("add", str(k))
#	OpCount.op("add", str(k))
#	OpCount.op("square", str(k))
#	OpCount.op("square", str(k))
#	if r2 == 0:
#		r1 = 1
#		r2 = 0
#	else:
#		r1 = r1/r2
#		r2 = 1
#	return [r1, r2]



#####################################

############## Lemma 3 ##############

#####################################

# l is the order of the kernel group
#k_prime and k_2_prime are defined the same as in the overleaf i.e if k or k_2 is even then it's divided by 2

# This function returns true if lemma 3 applies and false otherwise

def Condition_lemma3(k_prime,k_2,k_2_prime,l) :
  if k_2 == l-1:
    return True
  if k_2==(l-1)/2 :
    if l%4==3 :
      return True
    elif valuation_2(k_prime)==valuation_2(k_2_prime)+1 :
      return True
  return False

#####################################

########### Lemma 3 primes ##########

#####################################

# Returns a list of primes for which lemma 3 is verified
# Remark : this is not the complete list of such primes

def list_primes_lemma(n):
  list_of_fields = []
  for l in primes(100000):
      K = GF(l)
      element = K.2
      if element.multiplicative_order() == l-1: # We only look for primes for which 2 is a generator
        list_of_fields.append(l)
  return(list_of_fields)

#       Version where we get the full kernel polynomial and then evaluate

#*****************************************************************************

# S0 : list of the generators of Galois orbits (x-coordinates)
# T : Kernel group (x-coordinates)
# G : Kernel group generator (Point)
# A : Elliptic curve Montgomery coefficient (in K)

def get_S0_true(Gx,A,l,k_prime): # Returns S0 when lemma 3 applies

  #print("G: {} \n A: {}\n l:{}\n k_prime: {}".format(G, A, l, k_prime))

  S_0 = [Gx] #Initialize S_0 with our first generator of Galois orbit
  x = Gx
  idx = (l-1)/(2*k_prime) # number of Galois orbits
  if COUNT_ALL:
      OpCount.op("mult", str(k))
      OpCount.op("add", str(k))
      OpCount.op("div", str(k))
  for i in range (1,idx) :
    x = xDBL(x,A)
    S_0.append(x)
  #print("S_0:{}".format(S_0))
  return S_0

def get_S0_T_false(G,A,l,k_prime,K,E): # Returns S0 and T when lemma 3 does not apply

  S_0=[G[0]] #Initialize S_0 with our first generator of Galois orbit
  T=[]
  x0=G #used to compute next genrator
  frob=K.frobenius_endomorphism()
  x1=G #used to compute Galois orbits
  for i in range (k) : # generates Galois orbit from a generator
    T.append(x1[0])
    x1=E(frob(x1[0]), frob(x1[1]))
    OpCount.op("frob", str(K.degree()))
    OpCount.op("frob", str(K.degree()))
  while len(S_O) < (l-1)/2 :
    while x0[0] in T : #Find next generator
      x0=x0+G
      OpCount.op("add", str(k))
    x1=x0
    S_0.append(x0[0])
    for i in range (k) : # generates Galois orbit from a generator
        T.append(x1[0])
        x1=E(frob(x1[0]), frob(x1[1]))
        OpCount.op("frob", str(K.degree()))
        OpCount.op("frob", str(K.degree()))
  return T, S_0

#####################################

######### Get kernel points #########

#####################################



def get_kernel_polynomial_points (p,K,G,A,l, k) : # Return T using the functions as defined above
  E = EllipticCurve(K,[0,A,0,1,0])
  k_prime = copy(k)
  q = pow(p,k)
  lg = int(log(k,2))
  for i in range (abs(lg)) :
      if COUNT_ALL:
          OpCount.op("square", str(k))
  Gx=G[0]
  if k % 2 == 0:
    k_prime = k/2
    if COUNT_ALL:
        OpCount.op("div", str(k))
    K_even=GF((p,k_prime),'x')
    Gx=K_even(Gx)
    A= K_even(A)
  K2 = GF(l)
  element = K2(2)
  k_2 = element.multiplicative_order()
  k_2_prime = k_2
  if k_2 % 2 == 0:
    k_2_prime = k_2/2
    if COUNT_ALL:
        OpCount.op("div", str(k))
  #print("k_prime: {}\nk_2: {}\nk_2_prime: {}\nl: {}".format(k_prime,k_2,k_2_prime,l))
  if Condition_lemma3(k_prime,k_2,k_2_prime,l):
    #print("lemma yes")
    S_0 = get_S0_true(Gx,A,l,k_prime)
    T = copy(S_0)
    for i in range (1,k) :
    	for j in range (len(S_0)) :
            tmp = pow(p,i)
            for m in range (abs(int(log(i,2)))):
                if COUNT_ALL:
                    OpCount.op("square", str(k))
            T.append(pow(S_0[j],tmp))
            for m in range (abs(int(log(tmp,2)))):
                if COUNT_ALL:
                    OpCount.op("square", str(k))
  else :
  	#print("lemma no")
  	if k % 2 == 0:
  		T, S_0 = get_S0_T_false(G,A,l,k_prime,K_even,E)
  	else :
  		T, S_0 = get_S0_T_false(G,A,l,k_prime,K,E)

  print("S_O:",S_0)
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
  S_0=[G[0]]
  lam= get_frob_value()
  r=1
  T_L=[[0,0,0] for i in range (l-1)] #[b,m,n] where b is 1 if it has been encountered, m i is the Galois orbit generator, n is its frobenius powering
  for i in range (k) : # We check all points in Galois orbit generated by G
    power=pow(lam,i)%l
    for m in range (abs(int(log(i,2)))):
        if COUNT_ALL:
            OpCount.op("square", str(k))
    T_L[power-1]=[1,0,i]
  while len(S_0)< (l-1)/k : # While we haven't found all generators
    while T_L[r-1][0]==1 : # We look for the smallest multiple of G which is not in any Galois orbit previously calculated
      r+=1
      if COUNT_ALL:
          OpCount.op("add", str(k))
    previous=S_0[T_L[r-2][1]]*(pow(lam,T_L[r-2][2])%l)
    previous_2=S_0[T_L[r-3][1]]*(pow(lam,T_L[r-3][2])%l)
    if COUNT_ALL:
        OpCount.op("mult", str(k))
        OpCount.op("mult", str(k))
    new= xADD(previous,G,previous_2)
    S_0.append(new)
    for i in range (k) : # Check all points in the newly generated Galois orbit
      power=pow(lam,i)%l
      for m in range (abs(int(log(i,2)))):
          if COUNT_ALL:
              OpCount.op("square", str(k))
      T_L[(r*power)%l-1]=[1,r,i]
      if COUNT_ALL:
          OpCount.op("mult", str(k))
          OpCount.op("add", str(k))
    return(S_0)



def evaluate_from_G(p,k,G,A,l,P): # Returns the evaluation at P of the Kernel polynomial generated by G

    Gx=G[0]
    Px=P[0]
    k_prime = copy(k)
    if k % 2 == 0:
        k_prime = k/2
        if COUNT_ALL:
            OpCount.op("div", str(k))
        K_even = GF((p,k_prime),'x')
        Gx =K_even(Gx)
        Px=K_even(Px)
        A = K_even(A)
    K2 = GF(l)
    element = K2(2)
    k_2 = element.multiplicative_order()
    k_2_prime = k_2
    if k_2 % 2 == 0:
        k_2_prime = k_2/2
        if COUNT_ALL:
            OpCount.op("div", str(k))
    #print("k_prime: {}\nk_2: {}\nk_2_prime: {}\nl: {}".format(k_prime,k_2,k_2_prime,l))
    if Condition_lemma3(k_prime,k_2,k_2_prime,l): #Checks if lemma 3 applies or not
    	#print("lemma yes")
    	S_0=get_S0_true(Gx,A,l,k_prime)
    else :
    	#print("lemma no")
    	E = EllipticCurve(K,[0,K(A),0,1,0])
    	if k % 2 == 0:
    		S_0=get_S0_false(G,A,l,k_prime,K_even,E)
    	else :
    		S_0=get_S0_false(G,A,l,k_prime,K,E)
    print("S_O:",S_0)
    res=Px-S_0[0]
    if COUNT_ALL:
        OpCount.op("add", str(k))
    for i in range (1,len(S_0)) : #Multiplies all generators of Galois orbits
        res=res*(Px-S_0[i])
        if COUNT_ALL:
            OpCount.op("mult", str(k))
            OpCount.op("add", str(k))
    power=1
    for i in range (1,k) : # Frobenius powering
        power=power + pow(p,i)
        if COUNT_ALL:
            OpCount.op("add", str(k))
        for m in range (abs(int(log(i,2)))):
            if COUNT_ALL:
                OpCount.op("square", str(k))

    res=pow(res,power)

    for m in range (abs(int(log(power,2)))):
        if COUNT_ALL:
            OpCount.op("square", str(k))
    res= K(res)
    #print("so we have same ?",(P[0]-S_0[0])**(p**(k-1)),P[0]-S_0[0]**(p**(k-1)))
    return(res)


#       Get the isogeny

#****************************************************************************************

def compute_rationals_bis(E, T, l, kX, K, k, A): # Compute the rationnals function for getting isogeny from E with kernel points T

    a1,a2,a3, a4, a6= 0, A, 0, 1, 0 #coefficients
    S = range(1,(l+1)/2)
    if l > 1:
      assert max(S) == (l-1)/2


    YY = X^3+a2*X^2+a4*X+a6

    YYprime = YY.diff()

    Psi = prod(X-t for t in T) #kernel polynomial
    #Psi = kX(Psi) # in case l=1
    Psiprime = Psi.diff()
    Psiprimeprime = Psiprime.diff()
    xsum = sum(K(t) for t in T)
    Phi = 4*YY*(Psiprime^2-Psiprimeprime*Psi)-2*YYprime*Psiprime*Psi+(l*X-xsum)*Psi^2

    assert Phi == 4*(X^3+a2*X^2+a4*X+a6)*(Psiprime^2-Psiprimeprime*Psi)-2*(3*X^2+2*a2*X+a4)*Psiprime*Psi+(l*X-xsum)*Psi^2
    Phiprime = Phi.diff()
    Omega = Phiprime*Psi-2*Phi*Psiprime
    return Psi, Phi, Omega
