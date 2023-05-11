COUNT_ALL = 1
load("step1.sage")
load("Costello-Hisil.sage")
load("step2.sage")
from ctool import OpCount
from database_sql import Database_iso

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
    db_name = "CH.db"
    Database_iso.create_database(db_name)

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
            K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
            A = GF(p)(A)
            E = EllipticCurve(K, [0, A, 0, 1, 0])
            Q = E.random_point()
            P = point_finding(A,p,l,k)
            print('P is ', P, 'order P is ', order(P))
            xP = [P[0], P[2]]
            print("A, p, l, k, : {} {} {} {} ".format(A, p, l, k))
            # ------------------- COMPUTE COSTELLO-HISIL --------------------------------------
            OpCount.clean()
            kernel = kernel_points(xP, A, (l-1)/2)
            ch_iso = odd_iso(kernel,Q)
            str_iso_ = "Kernel: " + str(kernel) +" P: "+ str(Q) +  ": iso = " + str(ch_iso)
            dba = Database_iso(p, k, l, OpCount.field_op[str(k)], str_iso_, A, db_name)
            dba.insert()
            OpCount.clean()
