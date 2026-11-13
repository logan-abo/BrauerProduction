LoadPackage( "ctbllib" );

# 	Given a character table tbl, checks each character for the set of prime element orders distinct 
# from the underlying characteristic for which an element of that order in the underlying group 
# acts with a fixed point through the corresponding representation.

#FindAllowedFixedPoints := function(tbl)
#	local p, char, indices, dimFixeds, allowed, vector, allowedFixedPoints;
#	
#	allowedFixedPoints := [];
#
#	for char in Irr(tbl) do
#		vector := [];
#
#		for p in PrimeDivisors(Size(tbl)) do
#			if p = UnderlyingCharacteristic then continue; fi;
#
#			indices := Filtered([1..Size(char)], x->OrdersClassRepresentatives(tbl)[x] = p);
#			dimFixeds := List(indices, i->(1/p)*Sum([1..p], x->char[PowerMap(tbl, x, i)]));
#			allowed := Size(Filtered(dimFixeds, x->x>0))>0;
#
#			if allowed then Append(vector, [p]); fi;
#		od;
#		
#		Append(allowedFixedPoints, [vector]);
#	od;
#	
#	
#	return Set(allowedFixedPoints);
#end;


#Given an ordinary character table tbl, runs FindAllowedFixedPoints on that table and all the Brauer tables for the same group.

#Read("tables/2A13_mod11.g");
#2A13_mod11_table := mod11;
#
#Print(FindAllowedFixedPoints(2A13_mod11_table), "\n");
#
#precompute_table := BrauerTable("2.A13", 11);
#
#if precompute_table = fail then
#	Print("This table is not available in the ATLAS");
#
#else
#	Print(FindAllowedFixedPoints(precompute_table));
#fi;


LoadPackage("ctbllib");
LoadPackage("spinsym");

FindAllowedFixedPoints := function(tbl)
    local r, primes, char, p, indices, dimFixeds, vector, allowedFixedPoints;

    r := UnderlyingCharacteristic(tbl);
    primes := Filtered(PrimeDivisors(Size(tbl)), x -> x <> r);
    allowedFixedPoints := [];

    for char in Irr(tbl) do
        vector := [];
        for p in primes do
            indices := Filtered([1..Size(char)], x -> OrdersClassRepresentatives(tbl)[x] = p);
            dimFixeds := List(indices, i -> (1/p) * Sum([1..p], x -> char[PowerMap(tbl, x, i)]));
            if Size(Filtered(dimFixeds, x -> x > 0)) > 0 then
                Add(vector, p);
            fi;
        od;
        Add(allowedFixedPoints, vector);
    od;

    Print("=== ", Identifier(tbl), " (r = ", r, ") ===\n");
    Print("  Universal primes (Condition 1) : ", Intersection(allowedFixedPoints), "\n");
    for vector in Set(allowedFixedPoints) do
        Print("  A_chi = ", List(vector, q -> Concatenation(String(r), "-", String(q))),
              "  (missing: ", Difference(primes, vector), ")",
              "  <- chars: ", Filtered([1..Length(allowedFixedPoints)], i -> allowedFixedPoints[i] = vector), "\n");
    od;
    Print("\n");

    return Set(allowedFixedPoints);
end;

# 1. Verify on 2.A13 mod 11 (produced vs. SpinSym/ATLAS)
Read("tables/2A13_mod11.g");
FindAllowedFixedPoints(mod11);
FindAllowedFixedPoints(BrauerTable("2.A13", 11));

# 2. Run across all tables for 2.Alt(14)
ord14 := CharacterTable("2.Alt(14)");
FindAllowedFixedPoints(ord14);

Read("tables/2A14_mod11.g");
FindAllowedFixedPoints(mod11);

Read("tables/2A14_mod13.g");
FindAllowedFixedPoints(mod13);