// 
// Constructs Polymer example, for paper.
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


class Polymer {
public:
    int n;
    int d =3;
    int m;
    MatrixXi edges;
    VectorXd lengths; 

    VectorXd default_initial(int num, int dum=1){


        n = num;     // number of edges
        m = n+1;
        
        VectorXd x0;     // initial condition, a dn-dimensional vector
        double half_low;
        double height;
        double tmp;

        // Construct the edges
        edges.resize(m-2, 2);

        // polymer edges
        for (int i = 0; i < m-2; i++) {
            edges(i, 0) = i;
            edges(i, 1) = i + 1;
        }
        
        //cout << edges << endl;
        
        int pos = 0;
        
        
        // Construct the initial condition
        x0 = VectorXd :: Zero(n * d);

        // If n is odd, then the initial condition consists of (n-1)/2 congruent isosceles triangles between 0 and n/2
        // If n is even, then the initial consition consists of n/2 congruent isosceles triangles and a line of length 1 between 0 and n/2
        // e.g. If n = 3, then the coordinates of three points are (0.375, 0.927),(0.75, 0),(1.125, 0.927) 
        // If n = 4,  then the coordinates of four points are (0.25, 0.968), (0.5, 0), (0.75, 0.968), (1,0)
        if (n % 2 ==1){
            tmp = 2*n+2;
            half_low = n/tmp;
            height = sqrt(1-half_low*half_low);
            
        }else{
            tmp = 2*n;
            half_low = (n-2)/tmp;
            height = sqrt(1-half_low*half_low);
        }

        for(int i=0; i<n;i++){
            x0(i*d) = (i+1)*half_low;
            if (i % 2 ==0){
                x0(i*d+1) = height;
            }
        }
        
        lengths = VectorXd ::Zero(m);
        set_lengths_from_x(x0);

        return x0;


    }
    void set_param(const MatrixXi& edge, int num, int dim){
        n = num;
        m = num+1;
        d = dim;
        edges.resize(m-2,2);
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
        int total;
        total = n;
        for (int k=0; k < m-2; k++) {
            i = edges(k,0);
            j = edges(k,1);
            a = lengths(k);
            qvals(k) = -a*a;
            
            for (int l=0; l<d; l++) {
                d2 = x(i*d+l) - x(j*d+l);
                qvals(k) += d2*d2;
            }
        }
        qvals(m-2)= -lengths(m-2)*lengths(m-2);
        qvals(m-1) = -lengths(m-1)*lengths(m-1);
        // length(m-2) represents the constraints that length between the first point and (0,0) is 1
        // length(m-1) represents the constraints that length between the last point and (n/2,0) is 1
        for (int l=0; l<d; l++) {
            d2 = x(l);
            qvals(m-2) += d2*d2;
        }
        for (int l=0; l<d; l++) {
            if (l ==0 ){
                d2 = x((total-1)*d+l)- total/2;
            }else{
                d2 = x((total-1)*d+l);
            }
            qvals(m-1) += d2*d2;
        }
    }
    void dq(const VectorXd& x, SpMat& Dq) {    // jacobian of constraint
        std::vector<Trip> tripletList;  // holds row,col,val for constructing jacobian
        tripletList.reserve((m)*(2*d));
        int i,j;
        double val;
        int total;
        total = n;

        for (int k=0; k < m-2; k++) {
            i = edges(k,0);
            j = edges(k,1);
            for (int l=0; l<d; l++) {
                val = x(i*d+l) - x(j*d+l);
                tripletList.push_back(Trip(i*d+l,k,2*val));
                tripletList.push_back(Trip(j*d+l,k,-2*val));
            }
        }
        for (int l=0; l<d; l++) {
            val = x(0*d+l);
            tripletList.push_back(Trip(l,m-2,2*val));

        }
        for (int l=0; l<d; l++) {
            if (l==0){
                val = x((total-1)*d+l)-total/2;

            }else{
                val = x((total-1)*d+l);
            }
            tripletList.push_back(Trip((total-1)*d+l,m-1,2*val));

        }
        Dq.resize(d*total,m);
        Dq.setFromTriplets(tripletList.begin(), tripletList.end());
    }
    

};

