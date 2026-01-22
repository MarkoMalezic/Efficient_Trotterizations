# Trotter-Suzuki decompostion schemes data



The code implements the recursive formulae for high order Trotterizations, which allows the evaluation of their efficiency, if the scheme parameters are known.
If the parameters aren't known, the implementation includes a minimizer, that optimizes the theoretical error and provides efficient Trotter-Suzuki schemes.
In its current state the program supports schemes up to order 6.
The code is written in C++, with a dependency on the Eigen linear algebra library.
The building is done by CMake (see `CMakeLists.txt`) and the instructions for this can be found inside `build-run_instructions.txt`.
The class structure can be seen in the schematic below.
