

//#include <QCoreApplication>
#include <iostream>
#include <fstream>
#define EIGEN_USE_MKL_ALL
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
    if (argc<2){
        cout << "Usage: " << argv[0] << " " << " filename" << endl;
        return -1;
    }

    int n = 64*1024; // largest data length is 128k

    // read in a text file that contains a real matrix stored in column major format
    // but read it into row major format
    ifstream fin(argv[1]);
    if (!fin.is_open()){
        cout << "Failed to open " << argv[1] << endl;
        return -1;
    }

    VectorXcd data = VectorXcd::Zero(n); // define a complex vector

    // read fid into vector
    double cola=0, colb=0;
    for(int i = 0; i < n; i++){
        fin >> cola >> colb;
        data(i) = complex<double> (cola,colb);
        if (fin.eof())
            break;
    }

    // set parameters
    double f1 = double(0.2482); double f2 = double(0.2491);
    int p = 20; int M = (int) floor(64*1024/2); int L = 64*2*1024;
    VectorXd eta(11);
    for(int i=0; i<11;i++)
        eta(i) = (double) (0+i*5e-5);
    VectorXcd Xdata = VectorXcd::Zero(L);

    locape locape1;
    locape1.set_arg(data, f1, f2, p, M, L, eta);
    locape1.transf_data();
    clock_t time_start = 0, time_end = 0;
    time_start = clock();
    locape1.locape_calculate(argv[2]);
    time_end = clock();
    cout << " -> Elapsed time is "
       << (double)(time_end - time_start)/CLOCKS_PER_SEC<< " seconds." <<endl;

    // write spectrum to a file
    VectorXd real_spec = (locape1.XLoCapEV).real();
    VectorXd imag_spec = (locape1.XLoCapEV).imag();
    writedat(argv[2], &real_spec[0], &real_spec[locape1.fL.size()-1], &imag_spec[0], &imag_spec[locape1.fL.size()-1]);

    return 0;

}
