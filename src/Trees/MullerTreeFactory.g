#############################################################################
# MullerTreeFactory.g
#
# Müller-specific construction of the faithful weight-1 Brauer trees
# for the Schur double cover 2.A_n.
#############################################################################

Read("src/Trees/MullerTreeTemplate.g");



GenerateMullerTrees := function(n, p)
    local tbl, data, trees, d, filename, out;

    filename := Concatenation(
        "data/2A",
        String(n),
        "_mod",
        String(p),
        "_trees.g"
    );

    if IsExistingFile(filename) then
        Read(filename);

        if not IsBoundGlobal("trees") then
            Error(
                "tree data file did not define the variable trees: ",
                filename
            );
        fi;

        return ValueGlobal("trees");
    fi;

    tbl := CharacterTable(
        Concatenation(
            "2.Alt(",
            String(n),
            ")"
        )
    );;

    data := MullerWeightOneData(
        n,
        p,
        tbl
    );

    trees := List(
        data,
        d -> MullerResolveTree(tbl, d)
    );

    out := OutputTextFile(
        filename,
        false
    );

    SetPrintFormattingStatus(
        out,
        false
    );

    PrintTo(
        out,
        "trees := ",
        trees,
        ";\n"
    );

    CloseStream(out);

    return trees;
end;



#############################################################################
# Return the abacus data X_j(lambda).
#############################################################################

MullerAbacus := function(lambda, p)
    local X, part, x, y;

    X := List([1 .. p], j -> []);

    for part in lambda do
        x := QuoInt(part, p);
        y := RemInt(part, p);

        Add(X[y + 1], x);
    od;

    return X;
end;


#############################################################################
# Determine whether lambda is a p-bar core.
#############################################################################

IsPBarCore := function(lambda, p)
    local X, x, y;

    X := MullerAbacus(lambda, p);

    for y in [0 .. p - 1] do
        for x in X[y + 1] do

            if [x, y] <> [1, 0] then

                if not x - 1 in X[y + 1] then
                    return false;
                fi;

            fi;

        od;
    od;

    if 1 in X[1] then
        return false;
    fi;

    for y in [1 .. QuoInt(p - 1, 2)] do

        if 0 in X[y + 1] and 0 in X[p - y + 1] then
            return false;
        fi;

    od;

    return true;
end;


#############################################################################
# Return all strict p-bar cores of m.
#############################################################################

MullerPBarCores := function(m, p)
    local parts;

    if m = 0 then
        return [ [] ];
    fi;

    parts := Filtered(
        Partitions(m),
        lambda -> ForAll(
            [1 .. Length(lambda) - 1],
            i -> lambda[i] > lambda[i + 1]
        )
    );

    return Filtered(
        parts,
        lambda -> IsPBarCore(lambda, p)
    );
end;


#############################################################################
# Sign of a strict partition.
#############################################################################

MullerSign := function(lambda)

    if RemInt(
        Sum(lambda) - Length(lambda),
        2
    ) = 0 then
        return 1;
    fi;

    return -1;
end;


#############################################################################
# Construct lambda_i^+.
#############################################################################

MullerPlusPartition := function(lambda, p, y)
    local X, x, part, result, i;

    X := MullerAbacus(lambda, p);

    x := Maximum(X[y + 1]);

    part := x * p + y;

    result := ShallowCopy(lambda);

    i := Position(result, part);

    result[i] := part + p;

    Sort(
        result,
        function(a, b)
            return a > b;
        end
    );

    return result;
end;


#############################################################################
# Construct lambda_0.
#############################################################################

MullerZeroPartition := function(lambda, p)
    local result;

    result := ShallowCopy(lambda);

    Add(result, p);

    Sort(
        result,
        function(a, b)
            return a > b;
        end
    );

    return result;
end;


#############################################################################
# Construct lambda_i^-.
#############################################################################

MullerMinusPartition := function(lambda, p, y)
    local result;

    result := ShallowCopy(lambda);

    Add(result, y);
    Add(result, p - y);

    Sort(
        result,
        function(a, b)
            return a > b;
        end
    );

    return result;
end;


#############################################################################
# Convert a Müller partition into a SpinSym character position.
#############################################################################

MullerSpinPos := function(tbl, partition, split)
    local param, pos;

    if split = 0 then
        param := [2, partition];

    elif split = '+' then
        param := [2, [partition, '+']];

    elif split = '-' then
        param := [2, [partition, '-']];

    else
        Error("invalid split: ", split);
    fi;

    pos := Position(
        CharacterParameters(tbl),
        param
    );

    if pos = fail then
        Error(
            "could not find SpinSym character parameter ",
            param
        );
    fi;

    return pos;
end;


#############################################################################
# Resolve a Müller tree from partition labels to character positions.
#############################################################################

MullerResolveTree := function(tbl, data)
    local plus, zero, minus, x;

    plus := [];

    for x in data.plus do
        Add(
            plus,
            MullerSpinPos(
                tbl,
                x.partition,
                x.split
            )
        );
    od;

    zero := MullerSpinPos(
        tbl,
        data.zero.partition,
        data.zero.split
    );

    minus := [];

    for x in data.minus do
        Add(
            minus,
            MullerSpinPos(
                tbl,
                x.partition,
                x.split
            )
        );
    od;

    return MullerTreeTemplate(
        plus,
        zero,
        minus,
        data.epsilon,
        data.name
    );
end;


#############################################################################
# Construct Müller weight-1 block data.
#############################################################################

MullerWeightOneData := function(n, p, tbl)
    local cores, lambda,
          X, occupied,
          y, yList,
          plus, minus, zero,
          part, data,
          epsilon,
          results;

    if p mod 2 = 0 then
        Error("MullerWeightOneData requires odd p");
    fi;

    if n < p then
        Error("MullerWeightOneData requires n >= p");
    fi;

    cores := MullerPBarCores(n - p, p);

    results := [];

    for lambda in cores do

        X := MullerAbacus(lambda, p);

        occupied := [];

        for y in [1 .. p - 1] do

            if 0 in X[y + 1] then
                Add(occupied, y);
            fi;

        od;

        Sort(
            occupied,
            function(a, b)
                return
                    Maximum(X[a + 1]) * p + a
                    >
                    Maximum(X[b + 1]) * p + b;
            end
        );

        plus := [];

        for y in occupied do

            part := MullerPlusPartition(
                lambda,
                p,
                y
            );

            if MullerSign(part) = 1 then
                Add(
                    plus,
                    rec(
                        partition := part,
                        split := '+'
                    )
                );
            else
                Add(
                    plus,
                    rec(
                        partition := part,
                        split := 0
                    )
                );
            fi;

        od;

        zero := MullerZeroPartition(lambda);

        if MullerSign(zero) = 1 then
            zero := rec(
                partition := zero,
                split := '+'
            );
        else
            zero := rec(
                partition := zero,
                split := 0
            );
        fi;

        yList := [];

        for y in [1 .. QuoInt(p - 1, 2)] do

            if not (
                y in occupied
                or p - y in occupied
            ) then

                Add(yList, y);

            fi;

        od;

        Sort(
            yList,
            function(a, b)
                return a > b;
            end
        );

        minus := [];

        for y in yList do

            part := MullerMinusPartition(
                lambda,
                p,
                y
            );

            if MullerSign(part) = 1 then
                Add(
                    minus,
                    rec(
                        partition := part,
                        split := '+'
                    )
                );
            else
                Add(
                    minus,
                    rec(
                        partition := part,
                        split := 0
                    )
                );
            fi;

        od;

        epsilon := MullerSign(lambda);

        data := rec(
            name := String(lambda),
            plus := plus,
            zero := zero,
            minus := minus,
            epsilon := epsilon
        );

        Add(
            results,
            data
        );

    od;

    return results;
end;