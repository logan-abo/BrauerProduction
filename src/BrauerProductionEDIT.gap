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




BrauerTable_2AN_modP := function(
    n,
    p,
    ordtbl,
    quotienttbl,
    trees,
    fus
)
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
    
    
    allTreeChars := ValidateBrauerTreeData(trees, faithful, p);
    
    
    
    
    
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
ValidateBrauerTreeData := function(trees, faithful, p)
    
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





t := CharacterTable("2.Alt(14)");;
a := CharacterTable("A14");;



t := CharacterTable("2.Alt(14)");;
a := CharacterTable("A14");;


cpt := ClassParameters(t);;
cpa := ClassParameters(a);;

fus := List(cpt, x -> Position(cpa, [1, x[2]]));;

StoreFusion(t, fus, a);




mod13 := BrauerTable2A14Mod13(t, a, trees);;