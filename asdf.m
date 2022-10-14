Set_finding:=procedure(A,p,l,~list)

  K:=GF(p);
  E:= EllipticCurve([0,K!A,0,K!1,0]);
  j:=jInvariant(E);
  R<T>:=PolynomialRing(K);
  phi_l:=ClassicalModularPolynomial(l);
  phi_l_j:=Evaluate(phi_l,[j,T]);

  // res := false;
  // for poly in Factored_phi_l_j do 
  //   res:=Degree(poly[1]) eq 1 or res;
  // end for;
  // if res then
  if HasRoot(phi_l_j) then
    N:=#E;
    t:=p+1-N;
    p_l:=p mod l;
    t_l:= t mod l;
    M<x>:=PolynomialRing(Rationals());
    f:=x^2-t_l*x+p_l;
    K0:=GF(l);
    roots_f:=Roots(f, K0);

    printf "-- roots = %o\n", roots_f;
    for lambda in roots_f do
      k:=Order(K0!lambda[1]);
      printf "-- k is %o\n", k;
      K2:=GF(p^k);
      E1 := BaseExtend(E, K2);
      // E1:= EllipticCurve([0,K2!A,0,K2!1,0]);
      exponent := Exponent(AbelianGroup(E1));
      assert exponent mod l eq 0;
      rate:= exponent div l;
      repeat
        G:=rate*Random(E1);
      until Order(G) eq l and (Integers( )!lambda[1])*G eq E1![G[1]^p,G[2]^p];  // Check eigenvalue
      R2<T>:=PolynomialRing(K2);
      poly:=R2!1;
      point:=G;
      for i in [1 .. (l-1) div 2] do
        poly *:= (T-point[1]); 
        point := point +G;
      end for;
      printf "-- poly is %o\n", poly;
      poly:= ChangeRing(poly,K);
      printf "-- E is %o\n-- poly is %o\n", E, poly;
      E1, phi:=IsogenyFromKernel(E,poly); 
      subgroup := Kernel(phi);
      
      Append(~list, [* A,p,l,k,E,j,N,G,E1,phi,subgroup *]);
      end for;
    end if;
  
end procedure;      