load("ogvelu.sage")
load("step2.sage")
load("step1.sage")

# ------ Initial part of Tests --------




#------ SAGE CODE #-------
A = 33
p = 37
l = 2
k = 1

K = GF((p,k),'x')
F.<x> = GF(p**k)
kX.<X> = F[]
E = EllipticCurve(K,[0,A,0,1,0])

# test

print("test for point finding")
#P = point_finding(A, p, l, k)

#P = point_finding(A, p, l, k)
#print("test for optimized point finding")
N = get_h_k(A, p, k)
#print(N)

print("E.order(): {}".format(E.order()))
P = optimized_point_finding(A, p, k,l, N, K)
print("P:{} - order: {}".format(P, P.order()))
#print(E.cardinality()/5)
#print("Print P: {}".format(P))
#print(P.order())
#print("\t test for get_kernel_polynomial_points")
T = get_kernel_polynomial_points(p,K,P,A,l,k)
#print("res get_kernel_polynomial_points: {}".format(T))

#print("\t Start Rationals : compute_rationals_bis \t")
Psi, Phi, Omega = compute_rationals_bis(E, T, l, kX, k, A)
#print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))
#print("P: {} - order {}".format(P, P.order()))


a_inv= E.a_invariants()


#Psi, Phi, Omega = compute_rationals(E, P, kX, k, a_inv)
##print("Psi: {}\nPhi: {}\nOmega: {}".format(Psi, Phi, Omega))
E_phi, new_as = build_new_curve(kX, Psi, Phi, Omega, E, P, a_inv)
print("New EllipticCurve: {}".format(E_phi))
