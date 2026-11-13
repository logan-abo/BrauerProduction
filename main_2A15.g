Read("src/Utils/Arithmetic.g");
Read("src/Trees/BrauerTrees.g");
Read("src/Trees/MullerTreeFactory.g");
Read("src/Utils/Export.g");
Read("src/BrauerFactory.g");

mullerTrees15mod11 := GenerateMullerTrees(15, 11);
mod11 := SpinAltBrauerTable(15, 11);;

mullerTrees15mod13 := GenerateMullerTrees(15, 13);
mod13 := SpinAltBrauerTable(15, 13);;