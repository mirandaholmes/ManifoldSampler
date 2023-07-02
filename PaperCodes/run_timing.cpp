
// 
// Tests timing of matrix factorizations
//
//
// Created Dec 14, 2022, at 9pm PST, going about 1000km/hr above Iceland
// 
//
//

#include <iostream>
#include <iomanip>    // for setting precision
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Equations.hpp"
#include <random>
#include <math.h>
#include <ctime>             // for timing operations
#include <chrono>   // for high-precision timing
#include "Polygon.hpp"
#include "Lattice.hpp"
#include "Polymer.hpp"
#include "Matrix1.hpp"
using namespace Eigen;
using namespace std;
using namespace std::placeholders;  // for binding functions

#define PI  3.141592653589793238462643383279502884


typedef std::chrono::high_resolution_clock Clock;

// column ordering, for factoring sparse matrices
typedef AMDOrdering<int> AMD; 
typedef COLAMDOrdering<int> COLAMD;   // this one is the cheapest (but doesn't work with Cholesky decomposition)
typedef SparseLU<SpMat, COLAMD> LU;    // Solver for Newton step (faster for large n)
typedef SimplicialLLT<SpMat, Eigen::Lower, AMD> Chol;  // Cholesky decomposition; for computing tangent steps


// Main loop
int main(int argc,char *argv[]){   


    /* ------------------------------------------------------ */
    /*   Parse input                                          */
    /* ------------------------------------------------------ */

    int iseed = -1;   // only needed for polygon

    // parameters to set
    int nfac;        // number of points to use in factorization tests
    int example;      // which exmaple to run
    int n;

    // user didn't provide enough input arguments
    if(argc <= 4) {
        cout << "ERROR: must provide at least 4 input arguments" << endl;
        cout << "Usage: [program] [example] [seed-] [nfac] [n1]" << endl;
        cout << "     can optionally provide more n, in increasing order" << endl;
        return 0;
    } 

    // extract input arguments
    example = atoi(argv[1]);
    iseed = atoi(argv[2]);
    nfac = atoi(argv[3]);
    int numn = argc-4;   // number of n-values

    string datafile = "Output/timing_e" + to_string(example) + "_s" + to_string(iseed) + ".txt";


    cout << "Running example " << example << " with " << endl;
    cout << "  nfac = " << nfac << endl;
    cout << "  numn = " << numn << endl;    
    cout << "  seed = " << iseed << endl;
    cout << "  datafile = " << datafile << endl;



    /* ------------------------------------------------------ */
    /*   Loop through n-values                                */
    /* ------------------------------------------------------ */
    for(int i=0; i< numn; i++) {

        n = atoi(argv[4+i]);

        cout << "\n===== n = " << n << " =====" << endl;

        cout << " ------------- " << endl;
        cout << "Beginning sampling code" << endl;

        /* Set up Equations */
        Polymer poly1;
        Lattice poly2;
        Matrix1 poly3;
        Polygon poly4;
        Equations myeqns;

        // initalize equations
        VectorXd x0;
        if(example == 1) {
            x0=poly1.default_initial(n);
            myeqns.nvars = poly1.n*poly1.d;
            myeqns.neqns = poly1.m;
            myeqns.eval_q = std::bind(&Polymer::q, poly1,_1, _2);
            myeqns.eval_Dq = std::bind(&Polymer::dq, poly1,_1,_2);
        }
        if(example == 2) {
            x0=poly2.default_initial(n);
            myeqns.nvars = poly2.n*poly2.d;
            myeqns.neqns = poly2.m;
            myeqns.eval_q = std::bind(&Lattice::q, poly2,_1, _2);
            myeqns.eval_Dq = std::bind(&Lattice::dq, poly2,_1,_2);
        }
        if(example == 3) {
            x0=poly3.default_initial(n);
            myeqns.nvars = poly3.n*poly3.d;
            myeqns.neqns = poly3.m;
            myeqns.eval_q = std::bind(&Matrix1::q, poly3,_1, _2);
            myeqns.eval_Dq = std::bind(&Matrix1::dq, poly3,_1,_2);
        }
        if(example == 4) {
            x0=poly4.default_initial(n,iseed);
            myeqns.nvars = poly4.n*poly4.d;
            myeqns.neqns = poly4.m;
            myeqns.eval_q = std::bind(&Polygon::q, poly4,_1, _2);
            myeqns.eval_Dq = std::bind(&Polygon::dq, poly4,_1,_2);
        }


        cout << "nvars = " << myeqns.nvars << ", neqns = " << myeqns.neqns << endl;
        //cout << "x0 = " << x0.transpose() << endl;

        // Evaluate Jacobian
        SpMat Dq,R,A;
        myeqns.eval_Dq(x0,Dq);


        R = Dq;   
        

        // construct a random matrix with the same sparsity pattern
        // random number generators for coordinates
        mt19937  rng; ;
        int seed0 = random_device()(); // choose seed randomly based on device configurations
        rng.seed(seed0);    // initialize generator for random number stream
        normal_distribution<double>   normrnd; 

       // change values of Dq
        int ct = 0;
        for (int k=0; k<R.outerSize(); ++k) {
            for (SpMat::InnerIterator it(R,k); it; ++it) {
              it.valueRef() = normrnd(rng);
            }
        }

        A = R.transpose()*R;
        
        //cout << "Dq = \n" << R << endl;
        //cout << "R = \n" << R << endl;
        //cout << "A = \n" << A << endl;


        // set up matrix factorizations
        Chol chol;
        LU lu;

        chol.analyzePattern(A);
        lu.analyzePattern(A);

        
        // timing variables
        clock_t start_time1, start_time2;  
        double timeChol, timeLU, tsolveChol, tsolveLU;      


        start_time1 = clock();            // start timing now
        for(int i=0; i<nfac; i++) {
            chol.factorize(A);
        }
        timeChol = (clock() - start_time1 ) / (double) CLOCKS_PER_SEC;


        start_time2 = clock();
        for(int i=0; i<nfac; i++) {
            lu.factorize(A);
        }
        timeLU = (clock() - start_time2 ) / (double) CLOCKS_PER_SEC;
        
        cout << "time (Chol) = " << timeChol << endl;
        cout << "time (LU)   = " << timeLU << endl;

        
        // Now check equation solving time
        // create a random vector
        VectorXd b;
        b.resize(A.cols());
        for(int i=0; i<b.cols(); i++) {
            b(i) = normrnd(rng);
        }

        // timing of linear solve
        start_time1 = clock();            // start timing now
        for(int i=0; i<nfac; i++) {
            chol.solve(b);
        }
        tsolveChol = (clock() - start_time1 ) / (double) CLOCKS_PER_SEC;


        start_time2 = clock();            // start timing now
        for(int i=0; i<nfac; i++) {
            lu.solve(b);
        }
        tsolveLU = (clock() - start_time2 ) / (double) CLOCKS_PER_SEC;

        cout << "solve (Chol) = " << tsolveChol << endl;
        cout << "solve (LU) = " << tsolveLU << endl;

        
        /* Write diagnostics to file */
        fstream myfile;
        myfile.open(datafile, ios::app);
        myfile.precision(6);   // Set precision

        myfile << example << " ";  //1
        myfile << iseed << " ";  //2
        myfile << nfac << " ";  //3
        myfile << n << " ";   //4
        myfile << timeChol << " ";  //5
        myfile << timeLU << " ";  //6
        myfile << tsolveChol << " "; //7
        myfile << tsolveLU << " "; //8
        myfile << endl;
        myfile.close();
        


        /* Write diagnostics to screen */
        /*cout << "Done sampling! " << endl;
        cout << "  npts = " << npts << endl;
        cout << "  sig  = " << mysampler.sig << endl;
        cout << "  time (seconds) = " << mysampler.time << " seconds" << endl;
        cout << "  Avg acc        = " << mysampler.rej.acc_avg() << endl;
        cout << "  Avg niter      = " << mysampler.niteravg() << endl;
        cout << "     Avg niter (success) =  "  << mysampler.niteravg(1) << endl;
        cout << "     Avg niter (fail)    =  "  << mysampler.niteravg(-1) << endl;
        cout << "  Reject (project)      = " << mysampler.rej.project << endl;
        cout << "  Reject (metropolis)   = " << mysampler.rej.metropolis << endl;
        cout << "  Reject (reverse proj) = " << mysampler.rej.revproj << endl;
        cout << "  Reject (reverse val)  = " << mysampler.rej.revval << endl;
        cout << "  Reject (total)        = " << mysampler.rej.total() << endl;
        cout << "  # Nonzeros (Dq, Dq'*Dq, L): " << mysampler.sparsity().transpose() << endl;
        */

    }
}

