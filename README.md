# LoCapE_CPP
C++ implementation of LoCapE algorithm


## Implementation high-level discription

This is an C++ implementation of LoCapE algorithm, which is the localized Capon estimator {link}.  

### Prerequisites

* This implementation uses [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for matrix manipulation. 
* Intel Math Kernel Library ([MKL](https://software.intel.com/en-us/mkl) ) is used for mathematical computation. 
* OpenMP is used for parallelizing the codes.

### coding style

There are mainly three files: main.cpp, locape.cpp and locape.h
* main.cpp contains data input reading, pre-processing, excuting locape algorithm and writing output.
* locape.h defines a class (locape). It contains three functions. 
```
    void set_arg(VectorXcd &data_a, double f1_a, double f2_a, int p_a, int M_a, int L_a, VectorXd eta_a);
    void transf_data(); // transfer the original time-sequence data into frequency domain through predefined fourier matrix F
    void locape_calculate(); // calculate amplitude for given frequency and damping factor
```
* locape.cpp is the file that contains core implementation of the algorithm

## TODO
* Optimize the implementation with GPU

