// 
// Constructs Ngon example, for paper 
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


class Polygon {
public:
    int iseed;
    int n;
    int d=3;
    int extra;
    int m;
    MatrixXi edges;
    VectorXd lengths; 

    VectorXd default_initial(int num, int seed0=-1){
        n = num;
        iseed=seed0;
        extra = n;
        m = n+extra;
        VectorXd x0;
        edges.resize(m,2);

        // noise parameters
        double sigz = 0.5;  // noise on z-coordinates 
        double minr = 0.6;      // min and max radius (uniformly chosen)
        double maxr = 1;  

        for (int i = 0; i < n-1; i++) {
            edges(i, 0) = i;
            edges(i, 1) = i + 1;
        }
        edges(n-1, 0) = n-1;
        edges(n-1, 1) = 0;

        //random number generator
        VectorXd seed_list;
        seed_list.resize(5);
        seed_list << 45946313, 1087225662, 1635865672, 1791411095, 1407244209;
        int seed;

        if(iseed <0) {
            seed = random_device()(); // choose seed randomly based on device configurations
        } else {
          seed = seed_list(iseed);
        }
        cout << "seed1 = " << seed << endl;

        std::default_random_engine eng(seed);
        std::uniform_int_distribution<int> distr(0, n-1);


        int tmpindex1, tmpindex2;
        int curindex1, curindex2;
        int flag = 0;

        //this part is used to generate rest of the constraints (num = extra) 
        for (int i = 0; i < extra; i++) {    
            tmpindex1 = distr(eng);
            tmpindex2 = distr(eng);

            //check if the constraints is the edge
            if (abs(tmpindex1 - tmpindex2) == 1 || abs(tmpindex1-tmpindex2) == n-1 || tmpindex1 == tmpindex2) {
                i--;
                continue;
            }
            for (int j = 0; j < i; j++) {
                curindex1 = edges(n + j, 0);
                curindex2 = edges(n + j, 1);
                //check if the constraint is already existed
                if ((tmpindex1 == curindex1 && tmpindex2 == curindex2) || (tmpindex1 == curindex2 && tmpindex2 == curindex1)) {
                    flag = 1;
                    break;
                }
            }
            if (flag == 1) {    
                flag = 0;
                i--;
                continue;
            }
            edges(n + i, 0) = tmpindex1;
            edges(n + i, 1) = tmpindex2;
        }

        VectorXd seed_list2;
        seed_list2.resize(5);
        seed_list2 << 1673865141, 2078057879, 527807503, 1584050293, 495166393;


        // random number generators for coordinates
        mt19937  rng; ;
        if(iseed<0) {
            seed = random_device()(); // choose seed randomly based on device configurations
        } else {
          seed = seed_list2(iseed);
        }
        cout << "seed2 = " << seed << endl;

        rng.seed(seed);    // initialize generator for random number stream
        normal_distribution<double>   normrnd; 
        uniform_real_distribution<double> unifrnd;

        
        // Put points on a polygon
        x0 = VectorXd :: Zero(n*d);
        double inner_angle, outer_angle, angle_tmp;
        double tmp_angle = n;
        inner_angle = (tmp_angle-2)*180/tmp_angle;
        outer_angle = 180-inner_angle;
        for (int i = 1; i < n ; i++) {
            angle_tmp = 180 - i*outer_angle; 
            x0(i*d) = x0((i-1)*d) + cos(angle_tmp*PI/180);
            x0(i*d+1) = x0((i-1)*d+1) + sin(angle_tmp*PI/180);
            x0(i*d+2) = 0;
        }  

        // put centre of mass at origin
        double xc=0, yc=0;
        for (int i=1; i<n; i++) {
            xc += x0(i*d)/n;
            yc += x0(i*d+1)/n;
        }
        for (int i=1; i<n; i++) {
            x0(i*d) = x0(i*d) - xc;
            x0(i*d+1) = x0(i*d+1) - yc;
        }
        

        // Now add noise
        double r;
        for (int i = 1; i < n ; i++) {
            x0(i*d+2) = normrnd(rng)*sigz;   // random z-coordinate
            r = unifrnd(rng)*(maxr-minr) + minr;
            x0(i*d) = r*x0(i*d);
            x0(i*d+1) = r*x0(i*d+1);
        }

        // set lengths from current coordinates
        set_lengths_from_x(x0);

        //cout << "lengths "<< lengths.transpose() << endl;
        return x0;
    }


    void set_param(const MatrixXi& edge, int num, int dim, int ext ){
        n = num;
        extra = ext;
        m = extra+n;
        d = dim;
        edges.resize(m,2);
        edges = edge;
        lengths = VectorXd :: Zero(m);
    }
 
    void set_lengths_from_x(const VectorXd& x) {
        VectorXd tmp(m); 
        lengths = VectorXd::Zero(m);
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






