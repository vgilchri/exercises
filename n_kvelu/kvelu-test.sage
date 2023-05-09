COUNT_ALL = 1
load("step1.sage")
load("step2.sage")
from ctool import OpCount

# --------------- PARAMETERS ----------------------------------------------
A,p,l,k = 52, 131, 19, 1
N = 152


K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])
filename = "params.txt"

# ---------------------	RUN OUR ALGO ----------------------------------------



for line in open(filename, 'r'):
    line = line.split(',')
    k = int(line[0].replace("[", ""))
    A = int(line[1])
    p = int(line[2])
    l = int(line[3].replace("]", ""))

    if k%2 != 0 and k < 12 and l > 2:
        #print("A, p, l, k, : {} {} {} {} ".format(A, p, l, k))
        k_prime = copy(k)
        if k % 2 == 0:
            k_prime = k/2
        K2 = GF(l)
        element = K2(2)
        k_2 = element.multiplicative_order()
        k_2_prime = k_2
        if k_2 % 2 == 0:
             k_2_prime = k_2/2
        if Condition_lemma3(k_prime,k_2,k_2_prime,l): #TODO: Rerun when the condition lema is false
            print("A, p, l, k, : {} {} {} {} ".format(A, p, l, k))
            N = get_h_k(A,p,k)
            print("Ratio: {}".format(N))
            fact = factor(N)
            not_cofact = True
            for f in list(fact):
                if l in f:
                    if f[1] > 1:
                        not_cofact = False
            if not_cofact:
                K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
                A = GF(p)(A)
                E = EllipticCurve(K, [0, A, 0, 1, 0])
                Q = E.random_point()
                print("start finding point...")
                #G = point_finding(A,p,l,k)

                G = optimized_point_finding(A,p,k,l,N,K)
                print("finish finding point...")
                OpCount.clean()
                our_iso = evaluate_from_G(p,k,G,A,l,Q)
                print("A, p, l, k, Q, iso: {} {} {} {} {} {}".format(A, p, l, k, Q, our_iso))
