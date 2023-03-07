load("ogvelu.sage")
load("step2.sage")
load("step1.sage")

# ------ Initial part of Tests --------


K = GF((17,4),'x')
q = 17
k.<fieldgen> = GF(q)
kX.<X> = k[]
#------ SAGE CODE #-------
A= 13
l = 11
k = 5
p= 17
# test

print("test for point finding")
#P = point_finding(A, p, l, k)

K = GF((p,k),'x')
E = EllipticCurve(K,[0,A,0,1,0])
P = point_finding(A, p, l, k)
print(P)
print(P.order())
#print("test for optimized point finding")
N=get_h_k(A, p, k)
#print(N)
E = EllipticCurve(K,[0,A,0,1,0])


P = optimized_point_finding(A, p, k,l, N, K)
#print(E.cardinality()/5)
#print("Print P: {}".format(P))
#print(P.order())
print("\t test for get_kernel_polynomial_points")
T = get_kernel_polynomial_points(p,K,P,A,l,k)
print("res get_kernel_polynomial_points: {}".format(T))
print(len(T))
print("T[0]=",T[0])


print("test for elliptic curves constructions bis")
print("EllipticCurve: {}\nkX: {}\nP: {}".format(E, kX, P))
print("\t Start Rationals : compute_rationals_bis \t")
Psi, Phi, Omega = compute_rationals_bis(E, T, l, kX, k, A)
print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))
print("P: {} - order {}".format(P, P.order()))

tmp_a = K(A)
a4 = (1/(1-(tmp_a**2/3)))
a6 = tmp_a/(3*(2*tmp_a**2/9-1))

a_inv= E.a_invariants()


#Psi, Phi, Omega = compute_rationals(E, P, kX, k, a_inv)
##print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))
E_phi, new_as = build_new_curve(kX, Psi, Phi, Omega, E, P, a_inv)
print("New EllipticCurve: {}".format(E_phi))

isoge = E.isogeny(P)
print("Sage Iso: {}".format(isoge))
