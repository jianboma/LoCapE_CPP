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
## TODO
* Optimize the implementation with GPU

