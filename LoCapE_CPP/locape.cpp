
#include<locape.h>
#include <omp.h>
#define EIGEN_USE_MKL_ALL
#include "unsupported/Eigen/MatrixFunctions"
using namespace std;
using namespace Eigen;

void locape::set_arg(VectorXcd &data_a, double f1_a, double f2_a, int p_a, int M_a, int L_a, VectorXd eta_a)
{
    data = data_a;
    f1 = f1_a; f2 = f2_a;
    p = p_a; M = M_a; L = L_a;
    eta = eta_a;

    pn = int(floor(L/M)); // jump size
    num_col = 2*p+1;
    fv = Eigen::VectorXd::Zero(L);
    indcol = Eigen::VectorXi::Zero(num_col); // define a complex vector

    // assign fv
    double bin = (double) 1.0/L;
    #pragma omp parallel for
    for (int i=0;i<L;i++) {
        fv(i) = 0+i*bin;
    }
    // assign indcol
    #pragma omp parallel for
    for (int i=0;i<num_col;i++) {
        indcol(i) = -p*pn+i*pn;
    }

    // assign frequency range and fourier matrix
    jL = (int) floor(f1*L)+1;
    jH = (int) ceil(f2*L)+1;
    int mid;
    if (jH+indcol(indcol.size()-1)>L)
      mid = L-jL-indcol(0)+1;
    else {
        mid = jH+indcol(indcol.size()-1)-jL-indcol(0)+1;
    }
    indfreq = Eigen::VectorXi::Zero(mid);
    #pragma omp parallel for
    for (int i=0;i<mid;i++) {
            indfreq(i) = jL+indcol(0)+i-1;
    } // assign indfreq
    len_ind = (int) indfreq.size();
    N = (int) data.size(); //number of samples
    K = N-M+1; //number of smooth points

    // assign original data
    dataB = data.conjugate().colwise().reverse(); // conjugate and reverse slices
    XF = Eigen::VectorXcd::Zero(M); // define a complex vector to store foreward data
    XB = Eigen::VectorXcd::Zero(M); // define a complex vector to store backward data
    #pragma omp parallel for
    for (int i=0;i<M;i++) {
        XF(i) = data(i);
        XB(i) = dataB(i);
    } // assign first M elements to XF and XB

    // initialize YB, YF and F (fourier matrix)
    YF = MatrixXcd::Zero(len_ind, K);
    YB = MatrixXcd::Zero(len_ind, K);
    F = MatrixXcd::Zero(len_ind, M);
    #pragma omp parallel for
    for (int i1=0;i1<len_ind;i1++) {
        for (int i2=0;i2<M;i2++) {
            F(i1,i2) = complex<double> (cos(-2*pi*fv(indfreq(i1))*i2),sin(-2*pi*fv(indfreq(i1))*i2));
        }
    }
    // for output
    XLoCapE = Eigen::MatrixXcd::Zero(eta.size(),jH-jL+1);
    fL = Eigen::VectorXcd::Zero(jH-jL+1);
    #pragma omp parallel for
    for (int i=0;i<jH-jL+1;i++) {
        fL(i) = fv(i+jL);
    }
};



void locape::transf_data()
{
    YF.col(0) = F*XF;
    YB.col(0) = F*XB;
    // assign other elements based on previous elements to reduce computation
    Eigen::VectorXcd one_vec = Eigen::VectorXcd::Ones(len_ind);
    Eigen::VectorXcd vec1, vec2,vec3, vec4, vec5, vec6, vec7; complex<double> point1, point2, point3, point4;
//    #pragma omp parallel for
    for (int i=1;i<K;i++) {
        YF.col(i) = (YF.col(i-1)-data(i-1)*one_vec).cwiseProduct(F.col(1).conjugate())+data(M+i-1)*F.col(M-1);
        YB.col(i) = (YB.col(i-1)-dataB(i-1)*one_vec).cwiseProduct(F.col(1).conjugate())+dataB(M+i-1)*F.col(M-1);
    }
}

void locape::locape_calculate(){

    XLoCapEV = Eigen::VectorXcd::Zero(jH-jL+1);
    maxIndex =  Eigen::VectorXi::Zero(jH-jL+1);
    complex<double> mid_value = 0; // to record largest value in each damping factor iteration
    int mid_ind = 0; // to record index of largest value in each damping factor iteration
    // start with each frequency first, as the inverse of covariance can be used for each damping factor
    MatrixXcd X1 = Eigen::MatrixXcd::Zero(indcol.size(),YF.cols());
    MatrixXcd Xi1 = Eigen::MatrixXcd::Zero(indcol.size(),YF.cols());
    MatrixXcd Fsub = Eigen::MatrixXcd::Zero(indcol.size(),M);

    omp_set_dynamic(0);
    omp_set_nested(1);

//    #pragma omp parallel for
    for (int jl=jL;jl<jH;jl++) {
        // pick columns from transformed data
        for (int i1=0;i1<indcol.size();i1++) {
            X1.row(i1) = YF.row(jl+indcol(i1)-indfreq(0)-1);
            Xi1.row(i1) = YB.row(jl+indcol(i1)-indfreq(0)-1);
            Fsub.row(i1) = F.row(jl+indcol(i1)-indfreq(0)-1);
        }
        // calculate covariance matrix and its inverse
        MatrixXcd mat1 = X1.transpose().conjugate();
        MatrixXcd R = (X1*X1.transpose().conjugate()+Xi1*Xi1.transpose().conjugate())/2.0;
        MatrixXcd Rinv = R.inverse();
        mid_value = 0; // initialize
        mid_ind = 0; // initialize

        #pragma omp parallel for
        for (int neta=0;neta<eta.size();neta++) {
            // calculate scalars
            double Le = (1-exp(-2*eta(neta)*K))/(1-exp(-2*eta(neta)));
            double LeN = (1-exp(-2*eta(neta)*N))/(1-exp(-2*eta(neta)));
            if (eta(neta) == 0)
            {
                Le = K;
                LeN = N;
            }
            // calculate fourier vectors
            Eigen::VectorXcd a_vec1 = Eigen::VectorXcd::Zero(M);
//            #pragma omp parallel for
            for (int i=0;i<M;i++) {
                a_vec1(i) = complex<double> (exp(-1*eta(neta)*i)*cos((2*pi*fv(jl-1))*i),exp(-1*eta(neta)*i)*sin((2*pi*fv(jl-1))*i));
            }
            Eigen::VectorXcd a_vec = Fsub*a_vec1;

            // calculate least square vectors
            Eigen::VectorXcd s_vec = Eigen::VectorXcd::Zero(K);
//            #pragma omp parallel for
            for (int i=0;i<K;i++) {
                s_vec(i) = complex<double> (1.0/Le*exp(-1*eta(neta)*i)*cos((-2*pi*fv(jl-1))*i),1.0/Le*exp(-1*eta(neta)*i)*sin((-2*pi*fv(jl-1))*i));
            }

            // calculate the amplitude
            Eigen::VectorXcd mag = (a_vec.conjugate().transpose()*Rinv*X1*s_vec)/((a_vec.conjugate().transpose()*Rinv)*a_vec);
            XLoCapE(neta,jl-jL) = mag(0);
            if (abs(mag(0))> abs(mid_value)){
                mid_value = mag(0);
                mid_ind = neta;
            }
        }
        Eigen::VectorXcd vec1 = XLoCapE.col(jl-jL);
        XLoCapEV(jl-jL) = mid_value;
        maxIndex(jl-jL) = mid_ind;
    }
}

