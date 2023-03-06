#*****************************************************************************

#       Algorithm for Step 2 : Computing an \ell-Isogeny

#*****************************************************************************



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
      b=True
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
