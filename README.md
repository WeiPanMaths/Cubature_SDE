# Cubature SDE solver

Numerical work for the PhD thesis "A pipeline for high precision PDE computation" (https://ora.ox.ac.uk/objects/uuid:b9cc83db-94b5-4c67-bc3c-1a9b234854d6).

The code implements "cubature on Wiener space with recombination and nonlinear adaptive approximation" for pricing exotics. Current implementation uses basket options as test cases, and is benchmarked against QuantLib (https://www.quantlib.org/) implementation. For details on methods, see thesis.


# Install

Requirements are

- cmake
- Quantlib
- Intel MKL
- Boost

Before compilation, edit the top CMakeLists.txt lines 14 -- 21, so that correct paths to the required include/lib files are set.

To compile, do the following

    cd ${source_dir}
    mkdir build
    cmake ..
    cmake --build .

The current cmake build script generates Debug version only. The executable output pde01.exe should appear in "${source_dir}/build/Debug".

Compilation tested on Windows 10 using the following versions of the required packages
- cmake v3.20.3
- Quantlib v1.21
- Intel MKL 2019.5.281
- Boost v1_66_0
