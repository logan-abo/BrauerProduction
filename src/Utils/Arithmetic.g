# Utils

#############################################################################
# p-adic exponent of an integer
#############################################################################

PExponent := function(n, p)
    local e;

    e := 0;
    while RemInt(n, p) = 0 do
        n := QuoInt(n, p);
        e := e + 1;
    od;

    return e;
end;



SpinCharPos := function(tbl, c, split)
    local param, p;

    if split = 0 then
        # A single ordinary spin character:
        # GAP parameter = [2,c]
        param := [2,c];

    elif split = '+' then
        # First character of a split pair:
        # GAP parameter = [2,[c,'+']]
        param := [2,[c,'+']];

    elif split = '-' then
        # Second character of a split pair:
        # GAP parameter = [2,[c,'-']]
        param := [2,[c,'-']];

    else
        Error("split must be 0, '+', or '-'");
    fi;

    p := Position(CharacterParameters(tbl), param);

    if p = fail then
        Error("ordinary character parameter not found: ", param);
    fi;

    return p;
end;



#############################################################################
# Select a maximal independent collection of rows of a matrix.
#
# For a Brauer tree incidence matrix this gives l(B) independent rows.
#############################################################################

IndependentRows := function(M)
    local basic, rank, i, candidate, newrank, ncols;

    basic := [];
    rank := 0;
    ncols := Length(M[1]);

    for i in [1 .. Length(M)] do

        candidate := M{Concatenation(basic, [i])};
        newrank := RankMat(candidate);

        if newrank > rank then
            Add(basic, i);
            rank := newrank;
        fi;

        if rank = ncols then
            break;
        fi;
    od;

    if rank <> ncols then
        Error("decomposition matrix does not have full column rank");
    fi;

    return basic;
end;



#############################################################################
# Linear combination of a list of value vectors.
#
# coeffs = [a1,...,ar]
# rows   = [v1,...,vr]
#
# returns a1*v1 + ... + ar*vr
#############################################################################

LinearCombinationValueLists := function(coeffs, rows)
    local n, j, i;

    n := Length(rows[1]);

    return List(
        [1 .. n],
        j -> Sum(
            [1 .. Length(rows)],
            i -> coeffs[i] * rows[i][j]
        )
    );
end;
