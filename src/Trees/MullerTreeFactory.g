#############################################################################
# MullerTreeFactory.g
#
# Müller-specific construction of the faithful weight-1 Brauer trees
# for the Schur double cover 2.A_n.
#############################################################################

Read("src/Trees/MullerTreeTemplate.g");



GenerateMullerTrees := function(n, p)
    local tbl, data, trees, d, filename, out;

    if p mod 2 = 0 or not IsPrimeInt(p) then
        Error("GenerateMullerTrees requires odd prime p");
    fi;

    if n < p then
        Error("MullerWeightOneData requires n >= p");
    fi;

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

        trees := ValueGlobal("trees");

        if IsList(trees) then
            if Length(trees) > 0 then
                return trees;
            fi;

            if Length(MullerPBarCores(n - p, p)) = 0 then
                return trees;
            fi;
        fi;
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

            if x > 0 and [x, y] <> [1, 0] then

                if not (x - 1 in X[y + 1]) then
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

    parts := Filtered(
        parts,
        lambda -> IsPBarCore(lambda, p)
    );

    Sort(
        parts,
        function(a, b)
            local i;

            i := 1;

            while i <= Length(a) and i <= Length(b) do
                if a[i] <> b[i] then
                    return a[i] > b[i];
                fi;

                i := i + 1;
            od;

            return Length(a) > Length(b);
        end
    );

    return parts;
end;


#############################################################################
# Return a partition label without GAP's display padding.
#############################################################################

MullerPartitionName := function(partition)
    local name, i;

    if Length(partition) = 0 then
        return "[]";
    fi;

    name := Concatenation(
        "[",
        String(partition[1])
    );

    i := 2;

    while i <= Length(partition) do
        name := Concatenation(
            name,
            ",",
            String(partition[i])
        );

        i := i + 1;
    od;

    return Concatenation(name, "]");
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

    minus := [];

    if data.epsilon = 1 then

        if data.zero.split <> '+' then
            Error("MullerResolveTree: positive tree requires split zero");
        fi;

        for x in data.plus do
            if x.split <> 0 then
                Error("MullerResolveTree: positive tree branch must be unsplit");
            fi;

            Add(
                plus,
                MullerSpinPos(
                    tbl,
                    x.partition,
                    x.split
                )
            );
        od;

        zero := [
            MullerSpinPos(tbl, data.zero.partition, '+'),
            MullerSpinPos(tbl, data.zero.partition, '-')
        ];

        for x in data.minus do
            if x.split <> 0 then
                Error("MullerResolveTree: negative branch must be unsplit");
            fi;

            Add(
                minus,
                MullerSpinPos(
                    tbl,
                    x.partition,
                    x.split
                )
            );
        od;

    elif data.epsilon = -1 then

        for x in data.plus do
            if x.split <> '+' then
                Error("MullerResolveTree: positive branch must split");
            fi;

            Add(
                plus,
                [
                    MullerSpinPos(tbl, x.partition, '+'),
                    MullerSpinPos(tbl, x.partition, '-')
                ]
            );
        od;

        if data.zero.split <> 0 then
            Error("MullerResolveTree: negative tree requires unsplit zero");
        fi;

        zero := MullerSpinPos(
            tbl,
            data.zero.partition,
            data.zero.split
        );

        for x in data.minus do
            if x.split <> '+' then
                Error("MullerResolveTree: negative branch must split");
            fi;

            Add(
                minus,
                [
                    MullerSpinPos(tbl, x.partition, '+'),
                    MullerSpinPos(tbl, x.partition, '-')
                ]
            );
        od;

    else

        Error(
            "MullerResolveTree: epsilon must be +1 or -1"
        );

    fi;

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

    if p mod 2 = 0 or not IsPrimeInt(p) then
        Error("MullerWeightOneData requires odd prime p");
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

        zero := MullerZeroPartition(lambda, p);

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
            name := MullerPartitionName(lambda),
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