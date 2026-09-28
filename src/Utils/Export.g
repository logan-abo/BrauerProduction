# Export utils

############################################################################
##
## Stores character table as reloadable precomputed GAP table
##
############################################################################

StoreCharacterTable := function(filename, tbl, variableName)
    local out, i, info, j, class, component;

    out := OutputTextFile(filename, false);

    if out = fail then
        Error("could not open character table file: ", filename);
    fi;

    AppendTo(
        out,
        variableName,
        ":= function()\nlocal tbl, i;\ntbl:=rec();\n"
    );

    for i in [3, 6 .. Length(SupportedCharacterTableInfo)] do
        if Tester(SupportedCharacterTableInfo[i - 2])(tbl) then
            info := SupportedCharacterTableInfo[i - 2](tbl);

            if SupportedCharacterTableInfo[i - 1] = "Irr" then
                info := List(info, ValuesOfClassFunction);
            fi;

            if IsString(info) then
                info := ReplacedString(info, "\"", "\\\"");
                if '\n' in info then
                    info := SplitString(info, "\n");
                    AppendTo(out, "tbl.", SupportedCharacterTableInfo[i - 1], ":=\nConcatenation([\n");
                    for j in [1 .. Length(info) - 1] do
                        AppendTo(out, "\"", info[j], "\\n\",\n");
                    od;
                    AppendTo(out, "\"", Last(info), "\"\n]);\n");
                else
                    AppendTo(out, "tbl.", SupportedCharacterTableInfo[i - 1], ":=\n\"", info, "\";\n");
                fi;
            elif SupportedCharacterTableInfo[i - 1] = "ConjugacyClasses" then
                AppendTo(out, "tbl.ConjugacyClasses:=\n[\n");
                for class in info do
                    AppendTo(
                        out,
                        "ConjugacyClass(tbl.UnderlyingGroup,\n",
                        Representative(class),
                        "),\n"
                    );
                od;
                AppendTo(out, "];\n");
            else
                AppendTo(
                    out,
                    "tbl.",
                    SupportedCharacterTableInfo[i - 1],
                    ":=\n",
                    info,
                    ";\n"
                );
            fi;
        fi;
    od;

    if IsLibraryCharacterTableRep(tbl) then
        for component in SupportedLibraryTableComponents do
            if IsBound(tbl!.(component)) then
                info := tbl!.(component);
                AppendTo(out, "tbl.", component, ":=\n");
                if IsString(info) then
                    AppendTo(out, "\"", info, "\";\n");
                else
                    AppendTo(out, info, ";\n");
                fi;
            fi;
        od;
    fi;

    if HasConjugacyClasses(tbl) and HasSizesConjugacyClasses(tbl) then
        AppendTo(
            out,
            "for i in [1..Length(tbl.ConjugacyClasses)] do\n",
            "  SetSize(tbl.ConjugacyClasses[i],tbl.SizesConjugacyClasses[i]);\n",
            "od;\n"
        );
    fi;

    if IsLibraryCharacterTableRep(tbl) then
        AppendTo(out, "ConvertToLibraryCharacterTableNC(tbl);\n");
    else
        AppendTo(out, "ConvertToCharacterTableNC(tbl);\n");
    fi;

    AppendTo(
        out,
        "return tbl;\nend;\n",
        variableName,
        ":= ",
        variableName,
        "();\n"
    );

    CloseStream(out);
end;



############################################################################
##
## Prints Character Table as formatted text file
##
############################################################################


PrintTableText := function(filename, tbl, wrap)
    local out;

    if not IsBool(wrap) then
        Error("wrap must be true or false");
    fi;

    out := OutputTextFile(filename, false);

    SetPrintFormattingStatus(out, wrap);
    LogTo(out);
    Display(tbl);
    LogTo();

    CloseStream(out);

    Print("Printed character table as text in ", filename, "\n");
end;



#############################################################################
## 
## CSV helpers
##
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