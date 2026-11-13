LoadPackage("ctbllib");
LoadPackage("spinsym");

#############################################################################
# Normalize a table object (CharacterTable, record with .Irr, or matrix)
# into a list of irreducible character value lists.
#############################################################################
ExtractIrrMatrix := function(tblObj)
    if IsCharacterTable(tblObj) then
        return List(Irr(tblObj), ValuesOfClassFunction);
    elif IsRecord(tblObj) and IsBound(tblObj.Irr) then
        if IsList(tblObj.Irr) and Length(tblObj.Irr) > 0
           and IsClassFunction(tblObj.Irr[1]) then
            return List(tblObj.Irr, ValuesOfClassFunction);
        fi;
        return tblObj.Irr;
    elif IsMatrix(tblObj) or IsList(tblObj) then
        return tblObj;
    fi;
    Error("Unsupported table format for verification");
end;

#############################################################################
# Compare the produced 2.A13 mod p Brauer table against the standard GAP
# ATLAS library Brauer table (CharacterTable("2.A13") mod p).
#############################################################################
VerifyRedundantBrauerTable := function(prodTbl, p)
    local atlasOrd, atlasMod, prodMat, atlasMat,
          prodDegs, atlasDegs, transTbl, transMat, permProdMat;

    Print("=================================================================\n");
    Print("Verifying 2.A13 mod ", p, " against ATLAS library (CTblLib)\n");
    Print("=================================================================\n");

    atlasOrd := CharacterTable("2.A13");
    if atlasOrd = fail then
        Error("Could not load ordinary ATLAS table 2.A13 from CTblLib");
    fi;

    atlasMod := atlasOrd mod p;
    if atlasMod = fail then
        Error("Could not load Brauer table 2.A13 mod ", p, " from CTblLib");
    fi;

    prodMat  := ExtractIrrMatrix(prodTbl);
    atlasMat := ExtractIrrMatrix(atlasMod);

    # 1. Check dimensions (number of irreducible Brauer characters x p-regular classes)
    Print("Produced table dimensions : ",
          Length(prodMat), " x ", Length(prodMat[1]), "\n");
    Print("ATLAS library dimensions  : ",
          Length(atlasMat), " x ", Length(atlasMat[1]), "\n");

    if DimensionsMat(prodMat) <> DimensionsMat(atlasMat) then
        Print("[FAIL] Dimension mismatch for 2.A13 mod ", p, "!\n");
        return false;
    fi;
    Print("[PASS] Dimensions match.\n");

    # 2. Check character degree multisets (column 1 = identity class 1a)
    prodDegs  := Collected(List(prodMat,  row -> row[1]));
    atlasDegs := Collected(List(atlasMat, row -> row[1]));
    if prodDegs <> atlasDegs then
        Print("[FAIL] Character degree multisets differ for 2.A13 mod ", p, "!\n");
        Print("  Produced degrees : ", prodDegs, "\n");
        Print("  ATLAS degrees    : ", atlasDegs, "\n");
        return false;
    fi;
    Print("[PASS] Character degree multisets match.\n");

    # 3. Check centralizer orders and class representative orders (if CharacterTable)
    if IsCharacterTable(prodTbl) then
        if Collected(SizesCentralizers(prodTbl))
           <> Collected(SizesCentralizers(atlasMod)) then
            Print("[FAIL] Centralizer size multisets differ!\n");
            return false;
        fi;
        if Collected(OrdersClassRepresentatives(prodTbl))
           <> Collected(OrdersClassRepresentatives(atlasMod)) then
            Print("[FAIL] Class representative order multisets differ!\n");
            return false;
        fi;
        Print("[PASS] Centralizer sizes and class orders match.\n");

        # Full character table isomorphism check (Irr + power maps + centralizers)
        transTbl := TransformingPermutationsCharacterTables(prodTbl, atlasMod);
        if transTbl <> fail then
            Print("[PASS] Full Brauer CharacterTable isomorphism verified via ",
                  "TransformingPermutationsCharacterTables.\n");
            Print("  Row permutation    : ", transTbl.rows, "\n");
            Print("  Column permutation : ", transTbl.columns, "\n\n");
            return true;
        fi;
    fi;

    # 4. Matrix isomorphism check on Irr values up to row/column permutation
    # (handles class reordering and central involution lift swaps g <-> z*g)
    transMat := TransformingPermutations(prodMat, atlasMat);
    if transMat = fail then
        Print("[FAIL] Irr matrices are NOT equivalent under row/column permutation!\n");
        return false;
    fi;

    permProdMat := Permuted(
        List(prodMat, row -> Permuted(row, transMat.columns)),
        transMat.rows
    );
    if permProdMat <> atlasMat then
        Print("[FAIL] Transformed matrix did not equal ATLAS matrix!\n");
        return false;
    fi;

    Print("[PASS] Irr matrices are identical up to row/column permutation.\n");
    Print("  Row permutation    : ", transMat.rows, "\n");
    Print("  Column permutation : ", transMat.columns, "\n\n");
    return true;
end;

#############################################################################
# Load or compute the 2.A13 mod 11 and mod 13 tables, then run verification.
#############################################################################
if not IsExistingFile("tables/2A13_mod11.g")
   or not IsExistingFile("tables/2A13_mod13.g") then
    Read("redundant_computation.g");
else
    Read("tables/2A13_mod11.g");
    Read("tables/2A13_mod13.g");
fi;

verified11 := VerifyRedundantBrauerTable(mod11, 11);
verified13 := VerifyRedundantBrauerTable(mod13, 13);

if verified11 and verified13 then
    Print("OVERALL RESULT: Both 2.A13 mod 11 and 2.A13 mod 13 match the ATLAS library.\n");
else
    Error("OVERALL RESULT: Verification failed for one or more tables.");
fi;