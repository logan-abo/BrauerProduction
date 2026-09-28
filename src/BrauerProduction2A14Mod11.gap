LoadPackage("ctbllib");
LoadPackage("spinsym");

#############################################################################
#
#  Brauer tables of cyclic-defect blocks from Brauer trees.
#
#  Intended application:
#       2.A14 modulo 11
#
#############################################################################

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


#############################################################################
# Convert a Müller-style tree record into:
#
#   vertices      = lists of ordinary character positions at each vertex
#   edges         = pairs of vertex numbers
#   ordchars      = all ordinary character positions in this block
#
# The exceptional vertex is always vertex 1.
#
#############################################################################

TreeGraphData := function(tree)
    local vertices, edges, v, branch, previous, j,
          ordchars;

    vertices := [ tree.center ];
    edges := [];

    for branch in tree.branches do

        previous := 1;

        for j in [1 .. Length(branch)] do

            Add(vertices, [branch[j]]);
            v := Length(vertices);

            Add(edges, [previous, v]);

            previous := v;
        od;
    od;

    ordchars := Set(Concatenation(vertices));

    return rec(
        vertices := vertices,
        edges := edges,
        ordchars := ordchars
    );
end;


#############################################################################
# Construct the decomposition matrix from a cyclic Brauer tree.
#
# Rows correspond to ordinary irreducible characters.
# Columns correspond to Brauer irreducible characters.
#
#############################################################################

TreeDecompositionMatrix := function(tree)
    local g, D, v, e, vertex, pos, row,
          ordchars, rowOf, i;

    g := TreeGraphData(tree);

    ordchars := g.ordchars;

    D := NullMat(Length(ordchars), Length(g.edges));

    # Map ordinary-character position -> row number in D
    rowOf := rec();

    for i in [1 .. Length(ordchars)] do
        rowOf.(String(ordchars[i])) := i;
    od;

    for e in [1 .. Length(g.edges)] do

        for v in g.edges[e] do

            for pos in g.vertices[v] do

                D[rowOf.(String(pos))][e] := 1;
            od;
        od;
    od;

    return rec(
        matrix := D,
        graph := g
    );
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


#############################################################################
# Recover Brauer characters from one Brauer tree.
#
# ordtbl = ordinary table
# modtbl = p-regular Brauer-table skeleton
# tree   = Müller-style tree record
#
#############################################################################

BrauerCharactersFromTree := function(ordtbl, modtbl, tree)
    local td, D, ordchars, irr,
          rest, basic, subD, invD,
          phi, j, i, reconstructed;

    td := TreeDecompositionMatrix(tree);

    D := td.matrix;
    ordchars := td.graph.ordchars;

    irr := Irr(ordtbl);

    #
    # Restrict the ordinary irreducibles to p-regular classes.
    #
    rest := List(
        ordchars,
        i -> ValuesOfClassFunction(
            RestrictedClassFunction(irr[i], modtbl)
        )
    );

    #
    # Choose l(B) independent rows.
    #
    basic := IndependentRows(D);

    subD := D{basic};
    invD := Inverse(subD);

    #
    # Solve
    #
    #     restrictions = D * BrauerCharacters
    #
    phi := List(
        [1 .. Length(D[1])],
        j -> LinearCombinationValueLists(
            invD[j],
            rest{basic}
        )
    );

    #
    # Verify every row of the decomposition matrix.
    #
    for i in [1 .. Length(ordchars)] do

        reconstructed :=
            LinearCombinationValueLists(D[i], phi);

        if reconstructed <> rest[i] then
            Error(
                "Brauer tree reconstruction failed for ordinary character ",
                ordchars[i]
            );
        fi;
    od;

    return rec(
        brauerValues := phi,
        decompositionMatrix := D,
        ordinaryCharacters := ordchars,
        basicSet := basic
    );
end;


#############################################################################
# Construct the complete 2.A14 mod 11 table.
#
# ordtbl    = ordinary table of 2.A14
# quotienttbl = ordinary table of A14
# trees     = all faithful positive-defect 11-blocks
#
#############################################################################

BrauerTable2A14Mod11 := function(ordtbl, quotienttbl, trees)
    local p, modtbl, modquotient,
          ordirr, zpos,
          faithful, nonfaithful,
          allTreeChars, treeResults,
          defectZeroFaithful,
          resultValues,
          tr, trres, i,
          fusmod, qirr, phi,
          inflatedValues,
          chi;

    p := 11;

    #####################################################################
    # Check |2.A14| has 11-part exactly 11.
    #####################################################################

    if PExponent(Size(ordtbl), p) <> 1 then
        Error("expected 11 to occur to the first power in |2.A14|");
    fi;

    #####################################################################
    # Construct the 11-regular table heads.
    #####################################################################

    modtbl := CharacterTableRegular(ordtbl, p);

    if modtbl = fail then
        Error("could not construct the 11-regular table of 2.A14");
    fi;

    modquotient := BrauerTable(quotienttbl, p);

    if modquotient = fail then
        Error("could not construct the 11-Brauer table of A14");
    fi;

    ordirr := Irr(ordtbl);

    #####################################################################
    # Locate the central involution z.
    #####################################################################

    zpos := First(
        [1 .. NrConjugacyClasses(ordtbl)],
        i -> OrdersClassRepresentatives(ordtbl)[i] = 2
             and SizesConjugacyClasses(ordtbl)[i] = 1
    );

    if zpos = fail then
        Error("could not identify the central involution");
    fi;

    #####################################################################
    # Faithful ordinary characters:
    #
    #        chi(z) = -chi(1)
    #####################################################################

    faithful := Filtered(
        [1 .. Length(ordirr)],
        i -> ordirr[i][zpos] = -ordirr[i][1]
    );

    nonfaithful := Difference(
        [1 .. Length(ordirr)],
        faithful
    );

    #####################################################################
    # Validate the tree data.
    #####################################################################

    allTreeChars := [];

    for tr in trees do

        for i in TreeGraphData(tr).ordchars do

            if not i in faithful then
                Error(
                    "tree contains non-faithful character ",
                    i
                );
            fi;

            if PExponent(ordirr[i][1], p) <> 0 then
                Error(
                    "tree contains character ",
                    i,
                    " of non-positive defect"
                );
            fi;

        od;

        Append(
            allTreeChars,
            TreeGraphData(tr).ordchars
        );

    od;

    allTreeChars := Set(allTreeChars);

    #####################################################################
    # Every faithful character not in a positive-defect tree must have
    # defect zero.
    #####################################################################

    for i in faithful do

        if not i in allTreeChars then

            if PExponent(ordirr[i][1], p) <> 1 then
                Error(
                    "faithful character ",
                    i,
                    " is neither in a tree nor defect zero"
                );
            fi;

        fi;

    od;

    #####################################################################
    # Recover the faithful positive-defect Brauer characters from
    # the Müller Brauer trees.
    #####################################################################

    treeResults := [];

    for tr in trees do

        trres := BrauerCharactersFromTree(
            ordtbl,
            modtbl,
            tr
        );

        Add(treeResults, trres);

    od;

    #####################################################################
    # Start the final list of Brauer characters.
    #####################################################################

    resultValues := [];

    #####################################################################
    # NONFAITHFUL BRAUER CHARACTERS
    #
    # Since the kernel of
    #
    #             2.A14 -> A14
    #
    # has order 2, which is prime to 11, the nonfaithful 11-Brauer
    # characters of 2.A14 are precisely the inflations of the
    # irreducible 11-Brauer characters of A14.
    #
    # The ordinary fusion t -> a induces the fusion
    #
    #       modtbl -> modquotient.
    #####################################################################

    fusmod := GetFusionMap(modtbl, modquotient);

    if fusmod = fail then
        Error(
            "could not obtain induced fusion ",
            "2.A14 mod 11 -> A14 mod 11"
        );
    fi;

    qirr := Irr(modquotient);

    for phi in qirr do

        inflatedValues := List(
            fusmod,
            j -> phi[j]
        );

        Add(
            resultValues,
            inflatedValues
        );

    od;

    #####################################################################
    # FAITHFUL DEFECT-ZERO CHARACTERS
    #
    # A defect-zero ordinary character restricts to an irreducible
    # Brauer character on the p-regular classes.
    #####################################################################

    defectZeroFaithful := Filtered(
        faithful,
        i -> not i in allTreeChars
              and PExponent(ordirr[i][1], p) = 1
    );

    for i in defectZeroFaithful do

        chi := RestrictedClassFunction(
            ordirr[i],
            modtbl
        );

        Add(
            resultValues,
            ValuesOfClassFunction(chi)
        );

    od;

    #####################################################################
    # FAITHFUL POSITIVE-DEFECT CHARACTERS
    #
    # Append the irreducible Brauer characters obtained from each
    # cyclic-defect Brauer tree.
    #####################################################################

    for trres in treeResults do

        Append(
            resultValues,
            trres.brauerValues
        );

    od;

    #####################################################################
    # Sanity check: number of irreducible Brauer characters must equal
    # the number of p-regular conjugacy classes.
    #####################################################################

    if Length(resultValues) <> NrConjugacyClasses(modtbl) then

        Error(
            "wrong number of Brauer characters: got ",
            Length(resultValues),
            ", expected ",
            NrConjugacyClasses(modtbl)
        );

    fi;

    #####################################################################
    # Convert the value lists into GAP character objects.
    #####################################################################

    resultValues := List(
        resultValues,
        vals -> Character(modtbl, vals)
    );

    #####################################################################
    # Install the irreducible Brauer characters.
    #####################################################################

    SetIrr(modtbl, resultValues);

    return modtbl;
end;

t := CharacterTable("2.Alt(14)");;
a := CharacterTable("A14");;

A14SpinPos := function(tbl, c, split)
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



t := CharacterTable("2.Alt(14)");;

Tree_3 := rec(
    name := "[3]",
    center := [
        A14SpinPos(t, [11,3], '+'),
        A14SpinPos(t, [11,3], '-')
    ],
    branches := [
        [
            A14SpinPos(t, [10,3,1], 0),
            A14SpinPos(t, [9,3,2], 0),
            A14SpinPos(t, [7,4,3], 0),
            A14SpinPos(t, [6,5,3], 0)
        ],
        [
            A14SpinPos(t, [14], 0)
        ]
    ],
    exceptionalMultiplicity := 2
);

Tree_21 := rec(
    name := "[2,1]",
    center := [
        A14SpinPos(t, [11,2,1], 0)
    ],
    branches := [
        [
            A14SpinPos(t, [12,2], '+'),
            A14SpinPos(t, [13,1], '+')
        ],
        [
            A14SpinPos(t, [8,3,2,1], '+'),
            A14SpinPos(t, [7,4,2,1], '+'),
            A14SpinPos(t, [6,5,2,1], '+')
        ],
        [
            A14SpinPos(t, [12,2], '-'),
            A14SpinPos(t, [13,1], '-')
        ],
        [
            A14SpinPos(t, [8,3,2,1], '-'),
            A14SpinPos(t, [7,4,2,1], '-'),
            A14SpinPos(t, [6,5,2,1], '-')
        ]
    ],
    exceptionalMultiplicity := 1
);

a := CharacterTable("A14");;

cpt := ClassParameters(t);;
cpa := ClassParameters(a);;

fus := List(cpt, x -> Position(cpa, [1, x[2]]));;

StoreFusion(t, fus, a);

trees := [ Tree_3, Tree_21 ];

mod11 := BrauerTable2A14Mod11(t, a, trees);;


LogTo(
    "C:/Users/votlbob/OneDrive - Texas State University/Desktop/PrimeGraph/mod11_display.txt"
);;

Display(mod11);

LogTo();;


#############################################################################
# CSV helper
#############################################################################

CSVQuote := function(s)
    s := String(s);
    s := ReplacedString(s, "\"", "\"\"");
    return Concatenation("\"", s, "\"");
end;


CSVLine := function(entries)
    local result, i;

    if Length(entries) = 0 then
        return "";
    fi;

    result := String(entries[1]);

    for i in [2 .. Length(entries)] do
        result := Concatenation(
            result,
            ",",
            String(entries[i])
        );
    od;

    return result;
end;


#############################################################################
# Export 2.A14 mod 11 Brauer table.
#
# Columns:
#   Original SpinSym ClassParameters(t), restricted to 11-regular classes.
#
# Rows:
#   Nonfaithful:
#       inflated A14 Brauer characters
#
#   Faithful defect-zero:
#       corresponding ordinary 2.A14 SpinSym character
#
#   Faithful positive-defect:
#       Müller block + Brauer-tree edge + incident ordinary characters
#
#############################################################################

Export2A14Mod11CSV_labeled := function(tbl, quotienttbl, trees, filename)
    local out, p,
          classParams, charParams, irr,
          regularPos,
          qtbl, qirr, qcharParams,
          faithful, zpos, treeChars,
          defectZero,
          localFusion,
          tr, td, edge,
          endpoints, endpointLabels,
          row, header,
          i, j, value,
          treecount, dzcount, qcount,
          treeResult;

    p := 11;

    #########################################################################
    # Original SpinSym information
    #########################################################################

    classParams := ClassParameters(tbl);
    charParams := CharacterParameters(tbl);
    irr := Irr(tbl);

    #########################################################################
    # 11-regular classes in the ORIGINAL table.
    #########################################################################

    regularPos := Filtered(
        [1 .. NrConjugacyClasses(tbl)],
        i -> RemInt(
            OrdersClassRepresentatives(tbl)[i],
            p
        ) <> 0
    );

    #########################################################################
    # Complete A14 Brauer table.
    #########################################################################

    qtbl := BrauerTable(quotienttbl, p);
    qirr := Irr(qtbl);
    qcharParams := CharacterParameters(quotienttbl);

    #########################################################################
    # Central involution
    #########################################################################

    zpos := First(
        [1 .. NrConjugacyClasses(tbl)],
        i -> OrdersClassRepresentatives(tbl)[i] = 2
             and SizesConjugacyClasses(tbl)[i] = 1
    );

    if zpos = fail then
        Error("could not identify central involution");
    fi;

    #########################################################################
    # Faithful ordinary characters
    #########################################################################

    faithful := Filtered(
        [1 .. Length(irr)],
        i -> irr[i][zpos] = -irr[i][1]
    );

    #########################################################################
    # Characters occurring in the Müller trees
    #########################################################################

    treeChars := [];

    for tr in trees do
        Append(
            treeChars,
            TreeGraphData(tr).ordchars
        );
    od;

    treeChars := Set(treeChars);

    #########################################################################
    # Faithful defect-zero characters
    #########################################################################

    defectZero := Filtered(
        faithful,
        i -> not i in treeChars
             and PExponent(irr[i][1], p) = 1
    );

    #########################################################################
    # Induced fusion on the 11-regular classes.
    #
    # We stored the ordinary fusion
    #
    #       2.A14 -> A14
    #
    # so GAP can obtain the corresponding modular fusion.
    #########################################################################

    localFusion := GetFusionMap(
        CharacterTableRegular(tbl, p),
        qtbl
    );

    if localFusion = fail then
        Error(
            "could not obtain 2.A14 -> A14 modular fusion"
        );
    fi;

    #########################################################################
    # Open CSV
    #########################################################################

    out := OutputTextFile(filename, false);

    #########################################################################
    # HEADER
    #########################################################################

    header := [
        CSVQuote("BrauerCharacter"),
        CSVQuote("Source"),
        CSVQuote("OrdinaryCharacterReference")
    ];

    for i in regularPos do
        Add(
            header,
            CSVQuote(
                String(classParams[i])
            )
        );
    od;

    AppendTo(
        out,
        CSVLine(header),
        "\n"
    );

    #########################################################################
    # 1. NONFAITHFUL BRAUER CHARACTERS
    #########################################################################

    qcount := 0;

    for i in [1 .. Length(qirr)] do

        qcount := qcount + 1;

        row := [
            CSVQuote(
                Concatenation("NF_", String(qcount))
            ),
            CSVQuote("Nonfaithful"),
            CSVQuote(
                String(qcharParams[i])
            )
        ];

        for j in [1 .. Length(regularPos)] do

            value := qirr[i][localFusion[j]];

            Add(
                row,
                CSVQuote(String(value))
            );

        od;

        AppendTo(
            out,
            CSVLine(row),
            "\n"
        );

    od;

    #########################################################################
    # 2. FAITHFUL DEFECT-ZERO BRAUER CHARACTERS
    #########################################################################

    dzcount := 0;

    for i in defectZero do

        dzcount := dzcount + 1;

        row := [
            CSVQuote(
                Concatenation("DZ_", String(dzcount))
            ),
            CSVQuote("Faithful defect-zero"),
            CSVQuote(
                String(charParams[i])
            )
        ];

        for j in regularPos do

            value := irr[i][j];

            Add(
                row,
                CSVQuote(String(value))
            );

        od;

        AppendTo(
            out,
            CSVLine(row),
            "\n"
        );

    od;

    #########################################################################
    # 3. FAITHFUL POSITIVE-DEFECT BRAUER CHARACTERS
    #########################################################################

    treecount := 0;

    for tr in trees do

        td := TreeGraphData(tr);

        #
        # Compute the Brauer characters for this entire tree once.
        #
        treeResult := BrauerCharactersFromTree(
            tbl,
            CharacterTableRegular(tbl, p),
            tr
        );

        for edge in [1 .. Length(td.edges)] do

            treecount := treecount + 1;

            endpoints := td.edges[edge];

            #################################################################
            # Identify the ordinary characters incident to this edge.
            #################################################################

            endpointLabels := List(
                endpoints,
                v -> Concatenation(
                    List(
                        td.vertices[v],
                        x -> String(charParams[x])
                    )
                )
            );

            row := [
                CSVQuote(
                    Concatenation(
                        "TREE_",
                        String(treecount)
                    )
                ),

                CSVQuote(
                    Concatenation(
                        "Faithful positive-defect; block ",
                        tr.name,
                        "; edge ",
                        String(edge)
                    )
                ),

                CSVQuote(
                    Concatenation(
                        endpointLabels[1],
                        " <-> ",
                        endpointLabels[2]
                    )
                )
            ];

            #################################################################
            # Brauer values for this edge.
            #
            # brauerValues[edge] is indexed in the same order as the
            # CharacterTableRegular(tbl,11) classes.
            #################################################################

            for j in [1 .. Length(regularPos)] do

                value := treeResult.brauerValues[edge][j];

                Add(
                    row,
                    CSVQuote(String(value))
                );

            od;

            AppendTo(
                out,
                CSVLine(row),
                "\n"
            );

        od;

    od;

    #########################################################################
    # Close
    #########################################################################

    CloseStream(out);

    Print(
        "Exported ",
        qcount,
        " nonfaithful + ",
        dzcount,
        " faithful defect-zero + ",
        treecount,
        " faithful tree Brauer characters.\n"
    );

end;

#############################################################################
# Export a character table in essentially the same data order as Display(),
# but without character/class labels.
#############################################################################

CSVQuoteSimple := function(s)
    s := String(s);
    s := ReplacedString(s, "\"", "\"\"");
    return Concatenation("\"", s, "\"");
end;


CSVLineSimple := function(entries)
    local result, i;

    if Length(entries) = 0 then
        return "";
    fi;

    result := String(entries[1]);

    for i in [2 .. Length(entries)] do
        result := Concatenation(
            result,
            ",",
            String(entries[i])
        );
    od;

    return result;
end;


ExportDisplayDataCSV := function(tbl, filename)
    local out,
          classNames,
          centralizers,
          irr,
          powerMaps,
          primes,
          i, j,
          row,
          pm,
          value;

    out := OutputTextFile(filename, false);

    classNames := ClassNames(tbl);
    centralizers := OrdersCentralizers(tbl);
    irr := Irr(tbl);

    #########################################################################
    # 1. Centralizer orders
    #
    # Display() prints these in factorized form. Here we export the actual
    # integer centralizer orders.
    #########################################################################

    row := [];

    for i in [1 .. Length(centralizers)] do
        Add(row, CSVQuoteSimple(centralizers[i]));
    od;

    AppendTo(
        out,
        CSVLineSimple(row),
        "\n"
    );

    #########################################################################
    # 2. Class names
    #########################################################################

    row := [];

    for i in [1 .. Length(classNames)] do
        Add(row, CSVQuoteSimple(classNames[i]));
    od;

    AppendTo(
        out,
        CSVLineSimple(row),
        "\n"
    );

    #########################################################################
    # 3. Power maps
    #
    # Display() prints the available power maps as:
    #
    #     2P  ...
    #     3P  ...
    #     5P  ...
    #
    # We omit the "2P", "3P", etc. labels and write only the map itself.
    #########################################################################

    powerMaps := ComputedPowerMaps(tbl);

    for i in [2 .. Length(powerMaps)] do
        if IsBound(powerMaps[i]) then

            pm := powerMaps[i];

            row := [];

            for j in [1 .. Length(pm)] do
                value := pm[j];

                if IsInt(value) then
                    Add(
                        row,
                        CSVQuoteSimple(
                            classNames[value]
                        )
                    );
                else
                    Add(
                        row,
                        CSVQuoteSimple(value)
                    );
                fi;
            od;

            AppendTo(
                out,
                CSVLineSimple(row),
                "\n"
            );
        fi;
    od;

    #########################################################################
    # 4. Irreducible Brauer characters
    #
    # No X.1, X.2, ... labels are included.
    #########################################################################

    for i in [1 .. Length(irr)] do

        row := [];

        for j in [1 .. Length(irr[i])] do
            Add(
                row,
                CSVQuoteSimple(
                    irr[i][j]
                )
            );
        od;

        AppendTo(
            out,
            CSVLineSimple(row),
            "\n"
        );
    od;

    CloseStream(out);

    Print(
        "Exported Display data to ",
        filename,
        "\n"
    );
end;


#Export2A14Mod11CSV_labeled(
#    t,
#    a,
#    trees,
#    "2A14_mod11.csv"
#);


#############################################################################
# Helpers
#############################################################################

CSVLine := function(entries)
    local result, i;

    if Length(entries) = 0 then
        return "";
    fi;

    result := String(entries[1]);

    for i in [2 .. Length(entries)] do
        result := Concatenation(
            result,
            ",",
            String(entries[i])
        );
    od;

    return result;
end;


# Convert 1,2,... into A,B,...,Z,AA,AB,...
LetterName := function(n)
    local result, r;

    result := "";

    while n > 0 do
        r := RemInt(n - 1, 26);
        result := Concatenation(
            String(CharInt(65 + r)),
            result
        );
        n := QuoInt(n - 1, 26);
    od;

    return result;
end;


#############################################################################
# Compact character values
#
# Rational/integer values:
#     0  -> .
#
# Irrational values:
#     A, B, C, ...
#
# Complex conjugates / stars are represented using the same convention
# as GAP's Display:
#     /A
#     *A
#############################################################################

CompactValueMap := function(tbl)
    local irr, i, j, value,
          representatives, names,
          result, rep, name,
          found, conj, star, conjstar;

    irr := Irr(tbl);

    representatives := [];
    names := [];

    #########################################################################
    # Find one representative for each irrationality orbit.
    #########################################################################

    for i in [1 .. Length(irr)] do
        for j in [1 .. Length(irr[i])] do

            value := irr[i][j];

            if not IsRat(value) then

                found := false;

                for rep in representatives do
                    if value = rep
                       or value = ComplexConjugate(rep)
                       or value = StarCyc(rep)
                       or value = ComplexConjugate(StarCyc(rep))
                    then
                        found := true;
                        break;
                    fi;
                od;

                if not found then
                    Add(representatives, value);
                    Add(
                        names,
                        LetterName(Length(names) + 1)
                    );
                fi;
            fi;
        od;
    od;

    #########################################################################
    # Convert one value to its compact Display-style representation.
    #########################################################################

    result := function(value)
        local k, rep0, nm;

        if value = 0 then
            return ".";
        fi;

        if IsRat(value) then
            return String(value);
        fi;

        for k in [1 .. Length(representatives)] do

            rep0 := representatives[k];
            nm := names[k];

            if value = rep0 then
                return nm;
            fi;

            if value = ComplexConjugate(rep0) then
                return Concatenation("/", nm);
            fi;

            if value = StarCyc(rep0) then
                return Concatenation("*", nm);
            fi;

            if value = ComplexConjugate(StarCyc(rep0)) then
                return Concatenation("/*", nm);
            fi;
        od;

        # Should never be reached.
        return String(value);
    end;

    return result;
end;


#############################################################################
# Export the information displayed by Display(tbl), without labels.
#############################################################################

ExportDisplayDataCSV := function(tbl, filename)
    local out,
          centralizers,
          primes,
          factors,
          factorData,
          classNames,
          powerMaps,
          irr,
          compact,
          i, j, k,
          row,
          value,
          pm;

    out := OutputTextFile(filename, false);

    centralizers := SizesCentralizers(tbl);
    classNames := ClassNames(tbl);
    powerMaps := ComputedPowerMaps(tbl);
    irr := Irr(tbl);

    compact := CompactValueMap(tbl);

    #########################################################################
    # 1. Centralizer orders, factorized exactly in the style of Display()
    #
    # Example:
    #
    # 2,11,11,10,...
    # 3,5,5,4,...
    # 5,2,2,2,...
    #
    #########################################################################

    primes := Set(
        Concatenation(
            List(
                centralizers,
                x -> List(
                    Collected(FactorsInt(x)),
                    y -> y[1]
                )
            )
        )
    );

    Sort(primes);

    for i in primes do

        row := [];

        for j in centralizers do

            factors := Collected(FactorsInt(j));

            k := Position(
                List(factors, x -> x[1]),
                i
            );

            if k = fail then
                Add(row, ".");
            else
                Add(
                    row,
                    String(factors[k][2])
                );
            fi;
        od;

        AppendTo(
            out,
            CSVLine(row),
            "\n"
        );
    od;

    #########################################################################
    # 2. Class names
    #
    # Example:
    # 1a,2a,4a,2b,...
    #########################################################################

    row := [];

    for i in classNames do
        Add(row, i);
    od;

    AppendTo(
        out,
        CSVLine(row),
        "\n"
    );

    #########################################################################
    # 3. Power maps
    #
    # Example:
    # 1a,1a,2a,1a,...
    #
    # No "2P", "3P", etc. labels are written.
    #########################################################################

    for i in [1 .. Length(powerMaps)] do

        if IsBound(powerMaps[i]) then

            pm := powerMaps[i];

            row := [];

            for j in pm do
                Add(
                    row,
                    classNames[j]
                );
            od;

            AppendTo(
                out,
                CSVLine(row),
                "\n"
            );
        fi;
    od;

    #########################################################################
    # 4. Irreducible Brauer character values
    #
    # No X.1, X.2, ... labels.
    #########################################################################

    for i in [1 .. Length(irr)] do

        row := [];

        for j in [1 .. Length(irr[i])] do
            value := compact(irr[i][j]);
            Add(row, value);
        od;

        AppendTo(
            out,
            CSVLine(row),
            "\n"
        );
    od;

    CloseStream(out);

    Print(
        "Exported Display data to ",
        filename,
        "\n"
    );
end;


#ExportDisplayDataCSV(
#    mod11,
#    "C:/Users/votlbob/OneDrive - Texas State University/Desktop/PrimeGraph/2A14_mod11_raw.csv"
#);
