#ifndef LOCAPE_H
#define LOCAPE_H

#endif // LOCAPE_H

#define EIGEN_USE_MKL_ALL
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "unsupported/Eigen/FFT"
#include "unsupported/Eigen/MatrixFunctions"
#include "Eigen/Dense"
#include "Eigen/Core"
#include <iomanip>
#include <omp.h>

using namespace std;
using namespace Eigen;

class locape{


public:
    VectorXcd data;
    double f1; double f2;
    int p; int M; int L;
    VectorXd eta;
    Eigen::MatrixXcd XLoCapE; // define a matrix to record locape estimate
    Eigen::VectorXcd XLoCapEV; // define a vector to record locape estimate
    VectorXcd fL; // the frequencies analysed
    VectorXi maxIndex;

    const double pi = 3.141592653589793;
    int pn; // jump size
    int N; // data size
    int num_col;
    int jL; int jH;
    int len_ind;
    int K; // number of smooth points

    Eigen::VectorXd fv;
    Eigen::VectorXi indcol; // define a complex vector
    Eigen::VectorXi indfreq; // define a complex vector
    Eigen::VectorXcd dataB; // define a complex vector to store backward data
    Eigen::VectorXcd XF; // define a complex vector to store foreward data
    Eigen::VectorXcd XB; // define a complex vector to store backward data
    Eigen::MatrixXcd YF; // define a complex matrix to store transformed foreward data
    Eigen::MatrixXcd YB; // define a complex matrix to store transformed backward data
    Eigen::MatrixXcd F; // define a fourier matrix

    void set_arg(VectorXcd &data_a, double f1_a, double f2_a, int p_a, int M_a, int L_a, VectorXd eta_a);
    void transf_data(); // transfer the original time-sequence data into frequency domain through predefined fourier matrix F
    void locape_calculate(const char* filename); // calculate

private:
//    float pi = 3.1415926;
//    int pn; // jump size
//    int N; // data size
//    int num_col;
//    int jL; int jH;
//    int len_ind;
//    int K; // number of smooth points

//    Eigen::VectorXf fv;
//    Eigen::VectorXi indcol; // define a complex vector
//    Eigen::VectorXi indfreq; // define a complex vector
//    Eigen::VectorXcf dataB; // define a complex vector to store backward data
//    Eigen::VectorXcf XF; // define a complex vector to store foreward data
//    Eigen::VectorXcf XB; // define a complex vector to store backward data
//    Eigen::MatrixXcf YF; // define a complex matrix to store transformed foreward data
//    Eigen::MatrixXcf YB; // define a complex matrix to store transformed backward data
//    Eigen::MatrixXcf F; // define a fourier matrix





};
