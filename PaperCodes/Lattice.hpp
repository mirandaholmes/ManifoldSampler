// 
// Constructs Lattice example, for paper. 
//
// Created June 17, 2022
// 
//
//

#include <iostream>
#include <iomanip>    // for setting precision
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Equations.hpp"
#include "ManifoldSampler.hpp"
#include <random>
#include <math.h>


using namespace Eigen;
using namespace std;
using namespace std::placeholders;  // for binding functions

#define PI  3.141592653589793238462643383279502884


class Lattice {
public:
    int n;
    int d=2;
    int m;
    MatrixXi edges;
    VectorXd lengths; 
    int sl;
    VectorXd default_initial(int num, int dum=1){
        int sidelength = (int) round(sqrt(num));
        n = sidelength*sidelength;            //number of vertices
        m = sidelength*(sidelength-1)*2;      // number of edges
        VectorXd x0;     // initial condition, a dn-dimensional vector
        
        // Construct the edges
        edges.resize(m, 2);

        //construct the initial condition for edges.
        int pos = 0;
        for (int i = 0; i < sidelength; i++) {

            for (int j = 0; j < sidelength - 1; j++) {
                edges(pos, 0) = sidelength * i + j;
                edges(pos, 1) = sidelength * i + j + 1;
                pos++;
            }
            for (int z = 0; z < sidelength - 1; z++) {
                edges(pos, 0) = sidelength * z + i;
                edges(pos, 1) = sidelength * (z + 1) + i;
                pos++;
            }
        }

        
        
            // Construct the initial condition
        x0.resize(n * d);
        
        double cur_x = 0;
        double cur_y = 0;
        x0(0) = cur_x;
        x0(1) = cur_y;
        
        for (int i = 0; i < sidelength; i++) {


            for (int j = 0; j < sidelength; j++) {
                int index = i * sidelength + j;
                x0(2 * index) = j;
                x0(2 * index + 1) = sidelength -i -1;


            }

        }
        lengths = VectorXd :: Zero(m);
        set_lengths_from_x(x0);
        return x0;

    }

    void set_param(const MatrixXi& edge, int sidelength, int dim){
        sl = sidelength;
        n = sidelength*sidelength;
        m = sidelength*(sidelength-1)*2;
        d = dim;
        edges.resize(m,2);
        edges = edge;
        lengths = VectorXd :: Zero(m);
    }
 
    void set_lengths_from_x(const VectorXd& x) {
        VectorXd tmp(m); 
        q(x,tmp);
        lengths = tmp.cwiseSqrt();
    }

    void q(const VectorXd& x, VectorXd& qvals) {    // constraint
        int i,j;
        double a, d2;
        qvals.resize(m);
        for (int k=0; k < m; k++) {
            i = edges(k,0);
            j = edges(k,1);
            a = lengths(k);
            qvals(k) = -a*a;
            for (int l=0; l<d; l++) {
                d2 = x(i*d+l) - x(j*d+l);
                qvals(k) += d2*d2;
            }
        }
    }
    void dq(const VectorXd& x, SpMat& Dq) {    // jacobian of constraint
       std::vector<Trip> tripletList;  // holds row,col,val for constructing jacobian
        tripletList.reserve(m*(2*d));
        int i,j;
        double val;
        for (int k=0; k < m; k++) {
            i = edges(k,0);
            j = edges(k,1);
            for (int l=0; l<d; l++) {
                val = x(i*d+l) - x(j*d+l);
                tripletList.push_back(Trip(i*d+l,k,2*val));
                tripletList.push_back(Trip(j*d+l,k,-2*val));
            }
        }
        Dq.resize(d*n,m);
        Dq.setFromTriplets(tripletList.begin(), tripletList.end());
    }

};


//=================================================
//          Energy Function
//=================================================
double EnergyFunctions::lattice_corner(const VectorXd& x,  const double k, const int sidelength) {

    double U = 0;
    int n = sidelength * sidelength;

    for (int i = 1; i < sidelength; i++) {
        int left;
        int right;
        int index;
        double explen = sidelength * sqrt(2)/2;
        if ( sidelength ==3 )
            explen = 2;
        
        double tmp;
        if (i % 2 == 1) {
            for (int j = 0; j < sidelength; j=j+1) {
                index = i * sidelength + j;
                left = index - sidelength - 1;
                right = index - sidelength + 1;
                
                if (left >= 0 && left < n && left != 0) {
                    double length = sqrt((x(2 * index) - x(2 * left)) * (x(2*index) - x(2*left)) + (x(2*index+1) - x(2*left+1)) * (x(2*index+1) - x(2*left+1)));
                    //cout << " index = " << index << "left = " << left  << endl;
                    U = U + (length - explen) * (length - explen);
                }
                if (right >= 0 && right < n && j!=sidelength-1) {
                    double length = sqrt((x(2 * index) - x(2 * right)) * (x(2 * index) - x(2 * right)) + (x(2 * index + 1) - x(2 * right + 1)) * (x(2 * index + 1) - x(2 * right + 1)));
                    //cout << " index = " << index  << "right = " << right << endl;
                    U = U + (length - explen) * (length - explen);
                }
            }

        }
        else {
            for (int j = 0; j < sidelength; j = j + 1) {
                index = i * sidelength + j;
                left = index - sidelength - 1;
                right = index - sidelength + 1;
                
                if (left >= 0 && left < n && left!=0) {
                    double length = sqrt((x(2 * index) - x(2 * left)) * (x(2 * index) - x(2 * left)) + (x(2 * index + 1) - x(2 * left + 1)) * (x(2 * index + 1) - x(2 * left + 1)));
                    //cout << " index = " << index << "left = " << left  << endl;
                    U = U + (length - explen) * (length - explen);
                }
                if (right >= 0 && right < n && j != sidelength - 1) {
                    double length = sqrt((x(2 * index) - x(2 * right)) * (x(2 * index) - x(2 * right)) + (x(2 * index + 1) - x(2 * right + 1)) * (x(2 * index + 1) - x(2 * right + 1)));
                    //cout << " index = " << index << "right = " << right << endl;
                    U = U + (length - explen) * (length - explen);
                }
            }
        }

    }


    U = k * U;
    return U;
}
