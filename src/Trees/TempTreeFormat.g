
# Store Brauer tree data


t := CharacterTable("2.Alt(14)");;
a := CharacterTable("A14");;


dir := String(DirectoryCurrent());


Read("../src/Utils/Arithmetic.g");
Read("../archive/2A14_mod11_trees.g");





treesFile := OutputTextFile("../data/2A14_mod11_trees.g", false);
SetPrintFormattingStatus(treesFile, false);

PrintTo(
    treesFile,
    "trees := ",
    trees,
    ";\n"
);