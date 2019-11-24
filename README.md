# LoCapE_CPP
C++ implementation of LoCapE algorithm


## Implementation high-level discription

This is an C++ implementation of LoCapE algorithm, which is the localized Capon estimator {link}.  

### Prerequisites

* This implementation uses [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for matrix manipulation. 
* Intel Math Kernel Library ([MKL](https://software.intel.com/en-us/mkl) ) is used for mathematical computation. 
* OpenMP is used for parallelizing the codes.

### Coding style

There are mainly three files: main.cpp, locape.cpp and locape.h
* main.cpp contains data input reading, pre-processing, excuting locape algorithm and writing output.
* locape.h defines a class (locape). It contains three functions. 
```
    void set_arg(VectorXcd &data_a, double f1_a, double f2_a, int p_a, int M_a, int L_a, VectorXd eta_a);
    void transf_data(); // transfer the original time-sequence data into frequency domain through predefined fourier matrix F
    void locape_calculate(); // calculate amplitude for given frequency and damping factor
```
* locape.cpp is the file that contains core implementation of the algorithm


### Compile the implementation
Eigen library is needed. Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.You can download it and put it in the same directory as main.cpp or other places (you'll need to add the path to compiler).


The implementation uses MKL. MKL lib need to be provided. Suppose your MKL lib is installed in MKL_DIR (e.g. MKL_DIR=/opt/intel), the compiling code is
```
g++ -m64 -I$MKL_DIR/lib/intel64 main.cpp locape.h locape.cpp -L$MKL_DIR/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -o LoCapE
```


### Input arguments

Input arguments are:

1. input signal (input file path)
2. output file path
3. data length for analysis ('N' in C-codes)
4. localization region size ('p' in C-codes)
5. frequency grid number interpolation ('pn' in C-codes)
6. starting point of frequency range ('f1' in C-codes)
7. ending point of frequency range ('f1' in C-codes)

### Compare matlab implementation with C++ implementation

Once built, you can validate the implementation with matlab implementation. The matlab m-files are in calib folder. You can try either exponential signal or simulated NMR signals.

An example code uses built C++ implementation is as follow

```
./LoCapE ./data/Simulated_NMR.txt ./data/spectrum_simulated_NMR.txt 65536 20 2 0.2465 0.2532
```
## TODO
* Optimize the implementation with GPU
* amplitude refinement

