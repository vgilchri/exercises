load('our-algs.sage')

A,p,l,k = 52, 131, 19, 1



K = GF(p^k, 'x') # need this in order to have field consisten with point_finding()
A = GF(p)(A)
E = EllipticCurve(K, [0, A, 0, 1, 0])

G =  E.point([100,99,1])#point_finding(A,p,l,k)
G = [G[0], G[2]]
print("Point G: {}".format(G))
eval_points = [E.point([47,112,1])]#generate_points(E, l-1)
print("Points to eval: {}".format(eval_points))
images = algorithm_1(G, eval_points, A, l)
print("algorithm_1 images: {}".format(images))
images = algorithm_1_using_alg3(G, eval_points, A, l)
print("algorithm_1_using_alg3 images: {}".format(images))
images = algorithm_1_using_alg4(G, eval_points, A, l)
print("algorithm_1_using_alg4 images: {}".format(images))
