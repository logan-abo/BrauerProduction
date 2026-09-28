# Export utils

############################################################################
##
## Stores character table as reloadable precomputed GAP table
##
############################################################################

StoreCharacterTable := function(filename, tbl)
    local out;

    out := OutputTextFile(filename, false);

    LogTo(out);
    PrintCharacterTable(tbl);
    LogTo();

    CloseStream(out);

    Print("Stored character table in ", filename, "\n");
end;



############################################################################
##
## Prints Character Table as formatted text file
##
############################################################################


PrintTableText := function(filename, tbl)
    local out;

    out := OutputTextFile(filename, false);

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