

//#include <QCoreApplication>
#include <iostream>
#include <fstream>
#define EIGEN_USE_MKL_ALL // use intel MKL for mathematical computation
#include "unsupported/Eigen/FFT"
#include "Eigen/Dense"
#include <omp.h>
#include "Eigen/Core"
#include <iomanip>
#include <vector>
#include<locape.h>
#include<time.h>



template<class InputIterator, class InputIterator2>
void writedat(const char* filename,
              InputIterator xbegin, InputIterator xend,
              InputIterator2 ybegin, InputIterator2 yend,
              int xprecision=5, int yprecision=5)
{
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; xbegin != xend and ybegin != yend; ++xbegin, ++ybegin)
    f << std::setprecision(xprecision) << *xbegin << '\t'
      << std::setprecision(yprecision) << *ybegin << '\n';
  f.close();
}



using namespace Eigen;
using namespace std;
int main(int argc, char** argv){

//    QCoreApplication a(argc, argv);
    // check for an argument
    if (argc<5){
        cout << "Usage: " << argv[0] << " " << " Input_filepath Output_filepath Input_length start_freq end_freq" << endl;
        return -1;
    }

    int N = atoi(argv[3]); // largest data length is 128k

    // read in a text file that contains a real matrix stored in column major format
    // but read it into row major format
    ifstream fin(argv[1]);
    if (!fin.is_open()){
        cout << "Failed to open " << argv[1] << endl;
        return -1;
    }

    VectorXcd data = VectorXcd::Zero(N); // define a complex vector

    // read fid into vector
    double cola=0, colb=0;
    for(int i = 0; i < N; i++){
        fin >> cola >> colb;
        data(i) = complex<double> (cola,colb);
        if (fin.eof())
            break;
    }

    // set parameters
//    double f1 = double(0.2482); double f2 = double(0.24821);
    double f1 = atoi(argv[4]); double f2 = atoi(argv[5]);
    int p = atoi(argv[4]); int M = static_cast<int>(floor(N/2)) ; int L = N*2;
    VectorXd eta(11);
    for(int i=0; i<11;i++)
        eta(i) = static_cast<double> (0+i*5e-5);
    VectorXcd Xdata = VectorXcd::Zero(L);

    locape locape1;
    locape1.set_arg(data, f1, f2, p, M, L, eta);
    locape1.transf_data();
    clock_t time_start = 0, time_end = 0;
    time_start = clock();
    locape1.locape_calculate();
    time_end = clock();
    cout << " -> Elapsed time is "
       << static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC<< " seconds." <<endl;

    // write spectrum to a file
    VectorXd real_spec = (locape1.XLoCapEV).real();
    VectorXd imag_spec = (locape1.XLoCapEV).imag();
    writedat(argv[2], &real_spec[0], &real_spec[locape1.fL.size()-1], &imag_spec[0], &imag_spec[locape1.fL.size()-1]);

    return 0;

}
