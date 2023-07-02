
// 
// Runs Ngon example for paper, various n
//
//
// Modified Dec 7, 2022
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
#include "Polygon.hpp"
#include "Lattice.hpp"
#include "Polymer.hpp"
#include "Matrix1.hpp"
using namespace Eigen;
using namespace std;
using namespace std::placeholders;  // for binding functions

#define PI  3.141592653589793238462643383279502884


// Uncomment to choose an example. Be sure to set "example" in code
//typedef Polymer T;  // example = 1
//typedef Lattice T;  // example = 2
//typedef Matrix1 T;    // example = 3
typedef Polygon T;  // example = 4


// Main loop
int main(int argc,char *argv[]){   

    // Equations parameters
    int example = 4;   
    


    /* ------------------------------------------------------ */
    /*   Parse input                                          */
    /* ------------------------------------------------------ */

    // parameters to set
    double acc_target;    // target acceptance ratio
    double sig;           // sigma_high, for acc, or sigma for running simulation
    int projmethod;       // 1 = newton, 2 = symmetric
    int iseed = -1;      // only needed for polygon
    int n;


    // user didn't provide enough input arguments
    if(argc <= 5) {
        cout << "ERROR: must provide at least 5 input arguments" << endl;
        cout <<"Usage: [program] [acc_target] [sigma_high] [projmethod] [seed] [n1]" << endl;
        cout << "     can optionally provide more n, in increasing order" << endl;
        cout << "     set acc_target=0 for just sampling. In this case, sigma_high = sigma for sampling." << endl;
        return 0;
    } 

    // extract input arguments
    acc_target = atof(argv[1]);
    sig = atof(argv[2]);
    projmethod = atoi(argv[3]);
    iseed = atoi(argv[4]);
    int numn = argc-5;   // number of n-values

    string datafile = "Output/data_ngon_a" + to_string((int)round(acc_target*100)) + "_proj" + to_string(projmethod) + "_seed" + to_string(iseed) + ".txt";


    cout << "Running example " << example << " with " << endl;
    cout << "  acc_target = " << acc_target << endl;
    cout << "  sigma      = " << sig << endl;
    cout << "  projmethod = " << projmethod << endl;
    cout << "  seed       = " << iseed << endl;
    cout << "  numn = " << numn << endl;
    cout << "  datafile = " << datafile << endl;


    /* ------------------------------------------------------ */
    /*   Parameters; don't generally change these            */
    /* ------------------------------------------------------ */

    // Sampler parameters
    int npts = 1e6;        // how many points to generate
    int npts_acc = 1e5;    // number of points to estimate acceptance ratio

    double tol = 1e-5;     // tolerance for Newton's method   (optional)
    int maxIter = 100;     // max no. of iterations for Newton's method (optional)
    int dsave = 100;       // how often to save statistics (optional; default = 1)
    int ifdebug = 0;       // 1 = debug, 0 = nodebug
    int nburn = 1e4;       // number of burn-in points for acceptance ratio estimation



    /* ------------------------------------------------------ */
    /*   Loop through n-values                                */
    /* ------------------------------------------------------ */
    for(int i=0; i< numn; i++) {

        n = atoi(argv[5+i]);

        cout << "\n===== n = " << n << " =====" << endl;


        // Projection string
        string projstr;
        if(projmethod == 1)  projstr = "Newton";
        else projstr = "Symmetric";

        cout << " ------------- " << endl;
        cout << "Beginning sampling code" << endl;
        cout << "example = " << example << endl;
        cout << "n = " << n << endl;    
        cout << "projection method = " << projstr << endl;

        /* Set up Equations */
        T poly;
        Equations myeqns;
        VectorXd x0;
        x0.resize(n*poly.d);
        if(example == 1 || example == 2 || example == 3) {
            x0 = poly.default_initial(n);   // initializes example: initial condition, all parameters in class
        }
        else {
            x0 = poly.default_initial(n,iseed);   // initializes example: initial condition, all parameters in class
        }


        // initalize equations
        myeqns.nvars = poly.n*poly.d;
        myeqns.neqns = poly.m;
        myeqns.eval_q = std::bind(&T::q, poly,_1, _2);
        myeqns.eval_Dq = std::bind(&T::dq, poly,_1,_2);

        // Optionally bind to other functions
        if(example == 2) {
            myeqns.energy = std::bind(&EnergyFunctions::lattice_corner,_1,5,(int) round(sqrt(n)));
        }
        //myeqns.evalstats = std::bind(&StatisticsFunctions::nostats,_1);  // turn off statistics; don't save anything


        /* Create a sampler object */
        ManifoldSampler mysampler(myeqns,x0,sig,tol,maxIter);


        // Change projection method
        if(projmethod == 1) {
            mysampler.projmethod = ManifoldSampler::cProjNewton;
        } else {
            mysampler.projmethod = ManifoldSampler::cProjSym;
        }

        // Turn on debugging
        if(ifdebug==1) {
            mysampler.ifdebug = ManifoldSampler::cYes;
        }
        

        // Find sigma for desired acceptance ratio
        // sets "sig" within mysampler object
        // final x becomes new initial condition x0 after each run
        if(acc_target > 1e-6) {
            cout << "\nFinding sigma for target acceptance a = " << acc_target << " ..."<< endl;
            mysampler.sig_finder(acc_target, 1.5*mysampler.sig, npts_acc, nburn); // added factor of 1.5 for next time since each example is different
            sig = mysampler.sig;   
        }
        


        /* Sample! */
        cout << "\nSampling...." << endl;
        mysampler.sample(npts,dsave);

        
        /* Write diagnostics to file */
        fstream myfile;
        myfile.open(datafile, ios::app);

        myfile << example << " ";  //1
        myfile << projmethod << " ";
        myfile << npts << " ";
        myfile << n << " ";   //4
        myfile << mysampler.sig << " ";  //5
        myfile << mysampler.time << " ";  //6 
        myfile << mysampler.rej.acc_avg() << " ";  //7
        myfile << mysampler.niteravg() << " "; //8
        myfile << mysampler.niteravg(1) << " ";  //9
        myfile << mysampler.niteravg(-1) << " "; //10
        myfile << mysampler.rej.project << " ";  //11
        myfile << mysampler.rej.metropolis << " "; //12
        myfile << mysampler.rej.revproj << " "; //13
        myfile << mysampler.rej.revval << " ";  //14
        myfile << mysampler.rej.total() << " ";  //15
        myfile << mysampler.sparsity().transpose() << " ";  //16
        myfile << endl;
        myfile.close();


        /* Write diagnostics to screen */
        cout << "Done sampling! " << endl;
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

    }
}

