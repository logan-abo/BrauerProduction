Read("src/Utils/Arithmetic.g");
Read("src/Trees/BrauerTrees.g");
Read("src/Trees/MullerTreeFactory.g");
Read("src/Utils/Export.g");
Read("src/BrauerFactory.g");

mullerTrees16mod11 := GenerateMullerTrees(16, 11);
mod11 := SpinAltBrauerTable(16, 11);;

mullerTrees16mod13 := GenerateMullerTrees(16, 13);
mod13 := SpinAltBrauerTable(16, 13);;