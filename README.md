# tetraWild-Kokkos
Static libraries in folder lib and header files in folder include are from fTetrawild (https://github.com/wildmeshing/fTetWild).

The static libraries are compliled on Linux Centos. You may need to build fTetrawild on your machine to get those static libraries.

Besides GMP library (https://gmplib.org/) is needed. To install

To build the project, simply run

sudo yum install gmp 

./run_installlibs_cmake.sh

make

To run the code, simply

./main

The code runs fine wiht g++ but fail with kokkos nvcc_wrapper

Outputs are shown in the output.txt file.

