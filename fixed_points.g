LoadPackage("ctbllib");
LoadPackage("spinsym");


show_char_indices := false;


# Prime characteristics r in {2, 3, 5, 7}
for r in [2, 3, 5, 7] do
    FindAllowedFixedPoints(BrauerTable(CharacterTable("2.Alt(13)"), r), show_char_indices);
od;

# Prime characteristics r in {11, 13} from precomputed tables
Read("tables/2A13_mod11.g");
FindAllowedFixedPoints(mod11, show_char_indices);

Read("tables/2A13_mod13.g");
FindAllowedFixedPoints(mod13, show_char_indices);



FindAllowedFixedPoints := function(tbl, showIndices)
    local r, primes, char, q, classIndices, dimFixeds, vector, allowedFixedPoints, charIndices;

    r := UnderlyingCharacteristic(tbl);
    primes := Filtered(PrimeDivisors(Size(tbl)), x -> x <> r);
    allowedFixedPoints := [];

    for char in Irr(tbl) do
        vector := [];
        for q in primes do
            classIndices := Filtered([1..Size(char)], x -> OrdersClassRepresentatives(tbl)[x] = q);
            dimFixeds := List(classIndices, i -> (1/q) * Sum([1..q], x -> char[PowerMap(tbl, x, i)]));
            if Size(Filtered(dimFixeds, x -> x > 0)) > 0 then
                Add(vector, q);
            fi;
        od;
        Add(allowedFixedPoints, vector);
    od;

    Print(Identifier(tbl), "\n");

    for vector in Set(allowedFixedPoints) do
        charIndices := Filtered([1..Length(allowedFixedPoints)], i -> allowedFixedPoints[i] = vector);
        if showIndices then
            Print("  chi_", charIndices[1], " (all indices: ", charIndices, ") ==> fixed-point-free: ", Difference(primes, vector), "\n");
        else
            Print("  chi_", charIndices[1], " ==> fixed-point-free: ", Difference(primes, vector), "\n");
        fi;
    od;
    Print("\n");

    return Set(allowedFixedPoints);
end;