// 
// Example of how to use ManifoldSampler
//
// This example samples points on the surface of a d-dimensional ellipsoid. 
// All necessary equations are defined in this file. 
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


class Matrix1{
public:
    int n;
    int d=1;
    int m;
    MatrixXi edges;
    VectorXd lengths; 
    int sl;
    VectorXd default_initial(int num, int dum=1){
        int sidelength = (int) round(sqrt(num));
        n = sidelength*sidelength;            //number of vertices  
        m = (n+sidelength)/2;;      // number of edges
        VectorXd x0;     // initial condition, a dn-dimensional vector
        x0.resize(n *d );
        sl = sidelength;
        
        for (int i = 0; i < sidelength; i++) {
            for (int j = 0; j < sidelength; j++){
                if (i == j) {
                    x0(i * sidelength + j) = 1;
                }
                else {
                    x0(i * sidelength + j) = 0;
                }
            }
        }
        return x0;
    }

    // reshape v into a matrix. 
    // v must have length n
    MatrixXd mtx_of_vec(const VectorXd& v) {
        MatrixXd M;
        M.resize(sl,sl);
        for (int i = 0; i < sl; i++) {
            for (int j = 0; j < sl; j++){
                M(i,j) = v(i * sl + j);
            }
        }
        return M;
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
        int side = (int) round(sqrt(n));
        for (int j = 0; j < side; j++) {
            qvals(j) = -1;
            for (int i = 0; i < side; i++) {
                
                qvals(j) += x(i * side + j)*x(i*side+j);
                
            }
            
        }
        int cnt = side;
        for (int z = 0; z < side; z++) {
            for (int p = 0; p < z; p++) {
                qvals(cnt) = 0;
                for (int i = 0; i < side; i++) {
                    qvals(cnt) += x(i * side + z) * x(i * side + p);
                
                }
                cnt++;
            }
        }
    }
    void dq(const VectorXd& x, SpMat& Dq) {    // jacobian of constraint
        std::vector<Trip> tripletList;
        //tripleList.reserve()
        int side = (int) round(sqrt(n));

        for (int i = 0; i < side; i++) {
            for (int j = 0; j < side; j++) {
                tripletList.push_back(Trip(i * side + j, j, 2*x(i * side + j)));
                //cout << " first iter " << i * side + j << j << endl;
            }
            int column = side;
            for (int z = 0; z < side; z++) {
                for (int q = 0; q < z; q++) {
                    tripletList.push_back(Trip(i * side + z, column, x(i * side + q)));
                    tripletList.push_back(Trip(i * side + q, column, x(i * side + z)));
                    column++;
                    //cout << " second iter " << i * side + z << i * side + q << endl;
                
                }
            }
        }
        Dq.resize(n, (n  + sqrt(n)) / 2);
        Dq.setFromTriplets(tripletList.begin(), tripletList.end());

    }
};
