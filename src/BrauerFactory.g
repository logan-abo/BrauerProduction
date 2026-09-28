LoadPackage("ctbllib");
LoadPackage("spinsym");



#############################################################################
# Construct a complete modular character table for 2.A_n in characteristic p
# using:
#
#   - the Brauer table of A_n for the nonfaithful characters;
#   - defect-zero ordinary characters of 2.A_n;
#   - supplied Brauer trees for the faithful positive-defect blocks.
#
# INPUT:
#
#   n            = degree of the alternating group
#   p            = prime
#   ordtbl       = ordinary character table of 2.A_n
#   quotienttbl  = ordinary character table of A_n
#   trees        = Brauer trees for all faithful positive-defect blocks
#   fus          = ordinary class fusion map 2.A_n -> A_n
#
# ASSUMPTION:
#
#   p occurs to the first power in |2.A_n|.
#
#############################################################################
ConstructSpinAltBrauerTable := function(
    n,
    p,
    ordtbl,
    quotienttbl,
    trees
)
    local modtbl, modquotient,
          ordirr, zpos,
          faithful, nonfaithful,
          allTreeChars, treeResults,
          defectZeroFaithful,
          resultValues,
          tr, trres, i,
          fusmod, qirr, phi,
          inflatedValues,
          chi;

    CheckCyclicSubgroup(ordtbl, p);

    
    #########################################################################
    # Construct the p-regular table head of 2.A_n.
    #########################################################################

    modtbl := CharacterTableRegular(ordtbl, p);
    if modtbl = fail then

        Error(
            "could not construct the ",
            p,
            "-regular table of 2.A",
            n
        );

    fi;
    
    
    #########################################################################
    # Obtain the complete p-Brauer table of A_n.
    #
    # These give precisely the nonfaithful Brauer characters after inflation.
    #########################################################################
    
    modquotient := BrauerTable(quotienttbl, p);

    if modquotient = fail then

        Error(
            "could not construct the ",
            p,
            "-Brauer table of A",
            n
        );

    fi;


    ordirr := Irr(ordtbl);
    
    
    #########################################################################
    # Locate the central involution z.
    #
    # For the double cover 2.A_n, z is the unique class of order 2 and
    # class size 1.
    #########################################################################

    zpos := First(
        [1 .. NrConjugacyClasses(ordtbl)],
        i ->
            OrdersClassRepresentatives(ordtbl)[i] = 2
            and SizesConjugacyClasses(ordtbl)[i] = 1
    );

    if zpos = fail then

        Error(
            "could not identify the central involution of 2.A",
            n
        );

    fi;
    
    
    #########################################################################
    # Faithful ordinary characters:
    #
    #        chi(z) = -chi(1)
    #
    # Nonfaithful characters factor through A_n.
    #########################################################################
    
    faithful := Filtered(
        [1 .. Length(ordirr)],
        i -> ordirr[i][zpos] = -ordirr[i][1]
    );

    nonfaithful := Difference(
        [1 .. Length(ordirr)],
        faithful
    );
    
    allTreeChars := ValidateBrauerTreeData(
        trees,
        faithful,
        ordirr,
        p
    );

    for i in faithful do

        if not i in allTreeChars
           and PExponent(ordirr[i][1], p) <> 1
        then
            Error(
                "faithful character ",
                i,
                " is neither in a tree nor defect zero"
            );
        fi;

    od;

    treeResults := List(
        trees,
        tr -> BrauerCharactersFromTree(ordtbl, modtbl, tr)
    );

    resultValues := [];

    fusmod := GetFusionMap(modtbl, modquotient);

    if fusmod = fail then
        Error(
            "could not obtain induced fusion from 2.A",
            n,
            " to A",
            n,
            " in characteristic ",
            p
        );
    fi;

    qirr := Irr(modquotient);

    for phi in qirr do
        inflatedValues := List(fusmod, j -> phi[j]);
        Add(resultValues, inflatedValues);
    od;

    defectZeroFaithful := Filtered(
        faithful,
        i -> not i in allTreeChars
             and PExponent(ordirr[i][1], p) = 1
    );

    for i in defectZeroFaithful do
        chi := RestrictedClassFunction(ordirr[i], modtbl);
        Add(resultValues, ValuesOfClassFunction(chi));
    od;

    for trres in treeResults do
        Append(resultValues, trres.brauerValues);
    od;

    if Length(resultValues) <> NrConjugacyClasses(modtbl) then
        Error(
            "wrong number of Brauer characters: got ",
            Length(resultValues),
            ", expected ",
            NrConjugacyClasses(modtbl)
        );
    fi;

    resultValues := List(
        resultValues,
        vals -> Character(modtbl, vals)
    );

    SetIrr(modtbl, resultValues);

    return modtbl;
end;



#########################################################################
# Check that p occurs to the first power in |2.A_n|.
#
# This is required because the supplied Brauer-tree construction assumes
# cyclic defect groups of order p.
#########################################################################
CheckCyclicSubgroup := function(char_table, p)

    if PExponent(Size(char_table), p) <> 1 then

        Error(
            "p = ", p,
            " does not occur to the first power in provided group"
        );

    fi;
    
end;


#########################################################################
# Validate the supplied Brauer-tree data.
#########################################################################
ValidateBrauerTreeData := function(trees, faithful, ordirr, p)
    local allTreeChars, tr, i;
    
    allTreeChars := [];

    for tr in trees do

        for i in TreeGraphData(tr).ordchars do

            #################################################################
            # Every character appearing in a faithful Brauer tree must
            # actually be faithful.
            #################################################################

            if not i in faithful then

                Error(
                    "tree contains non-faithful character ",
                    i
                );

            fi;


            #################################################################
            # Every character in a tree must have positive p-defect.
            #
            # Since |G|_p = p, this means
            #
            #     p does not divide chi(1).
            #
            #################################################################

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
    
    return Set(allTreeChars);
    
end;


SpinAltBrauerTable := function(n, p)
    local tableFile, treeFile, variableName, trees,
          ordtbl, quotienttbl, classParams, quotientParams,
          fusion, modtbl;

    tableFile := Concatenation(
        "tables/2A",
        String(n),
        "_mod",
        String(p),
        ".g"
    );

    variableName := Concatenation("mod", String(p));

    if IsExistingFile(tableFile) then
        Read(tableFile);

        if not IsBoundGlobal(variableName) then
            Error(
                "precomputed table file did not define ",
                variableName,
                ": ",
                tableFile
            );
        fi;

        return ValueGlobal(variableName);
    fi;

    treeFile := Concatenation(
        "data/2A",
        String(n),
        "_mod",
        String(p),
        "_trees.g"
    );

    if not IsExistingFile(treeFile) then
        Error(
            "Brauer tree data file not found: ",
            treeFile,
            ". Generate this data before trying again."
        );
    fi;

    Read(treeFile);

    if not IsBoundGlobal("trees") then
        Error("tree data file did not define the variable trees: ", treeFile);
    fi;

    trees := ValueGlobal("trees");

    ordtbl := CharacterTable(
        Concatenation("2.Alt(", String(n), ")")
    );

    if ordtbl = fail then
        Error("could not retrieve the ordinary table of 2.Alt(", n, ")");
    fi;

    quotienttbl := CharacterTable(Concatenation("A", String(n)));

    if quotienttbl = fail then
        Error("could not retrieve the ordinary table of A", n);
    fi;

    classParams := ClassParameters(ordtbl);
    quotientParams := ClassParameters(quotienttbl);
    fusion := List(
        classParams,
        x -> Position(quotientParams, [1, x[2]])
    );

    if fail in fusion then
        Error("could not construct the ordinary class fusion 2.A", n, " -> A", n);
    fi;

    StoreFusion(ordtbl, fusion, quotienttbl);

    modtbl := ConstructSpinAltBrauerTable(
        n,
        p,
        ordtbl,
        quotienttbl,
        trees
    );

    StoreCharacterTable(tableFile, modtbl, variableName);

    return modtbl;
end;