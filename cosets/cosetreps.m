clear;

val2 := function(n)
    res := 0;
    while (n mod 2 eq 0) do 
        res := res + 1;
        n := Floor(n/2);
    end while;
    return res;
end function;

LemmaCheck := function(l,k)
    elt := GF(l)!2;
    k2 := Order(elt); // k_2
    k2prime := k2; // k_2_prime
    if k2 mod 2 eq 0 then 
        k2prime := Floor(k2/2);
    end if;
    k3 := k;
    if k mod 2 eq 0 then
        k3 := Floor(k/2); // k_prime
    end if;
    if k2 eq (l-1) then 
        return true;
    end if;
    if k2 eq (l-1)/2 then 
        if (l mod 4 eq 3) or (val2(k3) eq val2(k2prime)+1) then
            return true;
        end if;
    end if;
    return false;
end function;

function minimal_coset_reps(ell, k)
    // The sequence of minimal representatives for cosets of
    // the order-k subgroup of the multiplicative group
    assert ell ne 2 and IsPrime(ell);
    assert (ell - 1) mod k eq 0;
    // ell is not too big, so we can do this by brute force...
    // 1. find the subgroup of order k.
    R := Integers(ell);
    lambda := PrimitiveElement(R)^((ell - 1) div k);
    subgroup := [ lambda^i : i in [0..k-1] ];
    if R!-1 notin subgroup then
        subgroup cat:= [-e: e in subgroup];
    end if;
    // 2. find coset reps
    covered := {R|};
    cosetreps := [];
    cofactor := ((ell-1) div #subgroup);
    candidate := 0;
    while #cosetreps lt cofactor do
        candidate +:= 1;
        if candidate notin covered then
            representative := candidate;
            Append(~cosetreps, representative);
            covered join:= {candidate*g : g in subgroup};
        end if;
    end while;
    return cosetreps;
end function;

procedure print_data(ell)
    printf "Data for ell = %o\n", ell;
    print "    k : cosets : max(min reps) : ratio";
    for k in Divisors(ell - 1) do
        if k eq 1 then
            continue;
        end if;
        cosetreps := minimal_coset_reps(ell, k);
        N := Max(cosetreps);
        hope := (ell-1)/(IsOdd(k) select 2*k else k);
        ratio := RealField(5)!N/#cosetreps;
        printf "%5o : %6o : %-13o : %o\n", k, #cosetreps, N, ratio;
    end for;
end procedure;

procedure print_percents()
bad_l := 0;
total_l := 0;
for k in [1..12] do
for i in [3..10000] do 
    if IsPrime(i) then
        if (i-1) mod k eq 0 then
            if not(LemmaCheck(i,k)) then
                //i;
                total_l := total_l + 1;
                X := minimal_coset_reps(i,k);
                for i in [3..#X] do 
                    n := X[i] - X[i-1];
                    if not(n in X) then
                        //false;
                        bad_l := bad_l + 1;
                    end if;
                end for;
            end if;
        end if;
    end if;
end for;
"---------------------";
"k is ";
k;
"percent of bad l:";
Floor(bad_l / total_l * 100);
end for;
end procedure;

Lemma := function(l)
    if IsPrime(l) then
        if not(LemmaCheck(l,1)) then
            return true;
        end if;
    end if;
    return false;
end function;

procedure optimize_doubling(n)
// try to optimize doubling approach for bad l up to n
    for l in [3..n] do 
        // filter down to only bad l (can't use only doubling)
        if Lemma(l) then
            upper := Floor((l-1)/2);
            X := [2..upper]; // define the set S we want to compute
            root := GF(l)!2; // we start by doubling 2's
            elt := root;
            chain := []; // documents chain of computations
            R := [2]; // documents which roots we are doubling
            adds := 0;
            doubles := 0;
            while #X gt 0 do // we want to compute each value in X, so we remove them as we go, until nothing is left
                i := Max(Index(X,elt),Index(X,-1*elt)); // can be negative the elt too
                if i ne 0 then  // check whether the previous computation is still in the set we want
                    X := Remove(X,i);
                    chain := chain cat [elt];
                    elt := elt * 2; // double the previous elt
                    doubles := doubles + 1;
                else 
                    while not(root in X) do  // look for next necessary root
                        root := root + 1;
                    end while;
                    R := R cat [root];
                    elt := root;
                    adds := adds + 1;
                end if;
            end while;
            "for l = ",l," roots are ", R;
        end if;
    end for;
end procedure;
