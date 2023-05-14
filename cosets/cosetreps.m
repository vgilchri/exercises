function minimal_coset_reps(ell, k)
    // The sequence of minimal representatives for cosets of
    // the order-k subgroup of the multiplicative group
    assert ell ne 2 and IsPrime(ell);
    assert (ell - 1) mod k eq 0;
    // ell is not *too* big, so we can do this by brute force...
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
    for k in [1..11] do
    for i in [3..1000] do 
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
