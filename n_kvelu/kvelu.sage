load("ogvelu.sage")
load("step2.sage")
load("step1.sage")

# ------ Initial part of Tests --------

from ctool import OpCount


#------ SAGE CODE #-------
A = 8
p = 31
l = 19
k = 3

K = GF((p,k),'x')
K0=GF(p,'x')
F.<x> = GF(p**k)
kX.<X> = F[]

# test

print("test for point finding")
#P = point_finding(A, p, l, k)

E = EllipticCurve(K,[0,A,0,1,0])
E0 = EllipticCurve(K0,[0,A,0,1,0])
#P = point_finding(A, p, l, k)
#print("test for optimized point finding")
N=get_h_k(A, p, k)
#print(N)

print("E.order(): {}".format(E.order()))
G = optimized_point_finding(A, p, k,l, N, K)
print("G:{} - order: {}".format(G, G.order()))
#print(E.cardinality()/5)
#print("Print G: {}".format(G))
#print(G.order())

print("\t test for get_kernel_polynomial_points")
T = get_kernel_polynomial_points(p,K,G,A,l,k)
#print("res get_kernel_polynomial_points: {}".format(T))

print("\t Start Rationals : compute_rationals_bis \t")
Psi, Phi, Omega = compute_rationals_bis(E, T, l, kX, k, A)
#print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))
#print("G: {} - order {}".format(G, G.order()))


a_inv= E.a_invariants()


#Psi, Phi, Omega = compute_rationals(E, G, kX, k, a_inv)
##print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))
E_phi, new_as = build_new_curve(kX, Psi, Phi, Omega, E, G, a_inv)
print("New EllipticCurve: {}".format(E_phi))

OpCount.print_results()

#Find some prime st condition on lemma is wrong
print("\t test for lemma")
list_of_fields = []
for l1 in primes(100000):
    K = GF(l)
    element = K(2)
if (element.multiplicative_order() != l1-1) and (element.multiplicative_order() != l1-1/2):
    list_of_fields.append(l1)
print (list_of_fields)

print("\t test for evaluate_from_G")
P=E0.random_point()
print(P)
result = evaluate_from_G(p,k,G,A,l,P)
print(result)
result2=Psi(X=P[0])
print(result2)

