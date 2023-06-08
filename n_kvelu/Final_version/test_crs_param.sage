COUNT_ALL = 1
load("point_finding.sage")
load("Costello-Hisil.sage")
load("kernel_poly_evaluation.sage")
from ctool import OpCount
from database_sql import Database_iso

# --------------- PARAMETERS ----------------------------------------------
A,p,l,k = 52, 131, 19, 1
N = 152
p = 12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769
A = 10861338504649280383859950140772947007703646408372831934324660566888732797778932142488253565145603672591944602210571423767689240032829444439469242521864171


K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])
filename = "params.txt"

# ---------------------	RUN OUR ALGO ----------------------------------------



ll = [3, 5, 7, 11, 13, 17, 103, 523, 821]
for l in ll:
        k_prime = copy(k)
        if k % 2 == 0:
            k_prime = k/2
        K2 = GF(l)
        element = K2(2)
        a_l = element.multiplicative_order()
        a_l_prime = a_l
        if a_l % 2 == 0:
             a_l_prime = a_l/2
        K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
        A = GF(p)(A)
        E = EllipticCurve(K, [0, A, 0, 1, 0])
        Q = E.random_point()
        P = point_finding(A,p,k,l,K)
        OpCount.clean()
        #---- Common Costs -----
        our_iso = evaluate_from_G(p,k,P,A,l,Q)
        str_iso_ = "G: " + str(P) +" P: "+ str(Q) +  ": iso = " + str(our_iso)
        OpCount.print_results()
        print("A, p, l, k, Q, iso: {} {} {} {} {} {}".format(A, p, l, k, Q, our_iso))
        OpCount.clean()
        # ------------------- COMPUTE COSTELLO-HISIL --------------------------------------
        xP = [P[0], P[2]]
        OpCount.clean()
        kernel = kernel_points(xP, A, (l-1)/2)
        ch_iso = odd_iso(kernel,Q)
        str_iso_ = "Kernel: " + str(kernel) +" P: "+ str(Q) +  ": iso = " + str(ch_iso)
        OpCount.print_results()
        OpCount.clean()
        print("A, p, l, k, Q, ch_iso: {} {} {} {} {} {}".format(A, p, l, k, Q, ch_iso))
