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