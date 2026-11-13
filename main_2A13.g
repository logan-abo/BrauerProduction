Read("src/Utils/Arithmetic.g");
Read("src/Trees/BrauerTrees.g");
Read("src/Trees/MullerTreeFactory.g");
Read("src/Utils/Export.g");
Read("src/BrauerFactory.g");

mullerTrees13mod11 := GenerateMullerTrees(13, 11);
mod11 := SpinAltBrauerTable(13, 11);;

mullerTrees13mod13 := GenerateMullerTrees(13, 13);
mod13 := SpinAltBrauerTable(13, 13);;