
// 
// Runs each of the four examples in paper, with parameters set by user. 
//
// Change example by choosing a "typedef" then setting the example number correctly. 
// Then set parameters in code. 
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
//typedef Matrix1 T;  // example = 3
typedef Polygon T;    // example = 4


// Main loop
int main(int argc,char *argv[]){   

    /* ------------------------------------------------------ */
    /*   These are parameters you might want to change        */
    /* ------------------------------------------------------ */

    // Equations parameters
    int example = 4;    // make sure this is consistent with "typedef T"
    int n = 12;
    int iseed = 1;
    

    // Sampler parameters
    int npts = 2e2;        // how many points to generate
    double sig = 0.1;      // step size for sampler
    int projmethod = 2;    // Projection method: 1 = newton, 2 = symmetric (best)
    double tol = 1e-5;     // tolerance for Newton's method   (optional)
    int maxIter = 100;     // max no. of iterations for Newton's method (optional)
    int dsave = 2;         // how often to save statistics (optional; default = 1)
    int ifdebug = 0;       // 1 = debug, 0 = nodebug
    string datafile = "test.txt";       // file to which to save data
    //datafile = "";                    // set to "" to not save

    // Parameters for sigma-finder
    double acc_target = 0;  // target acceptance ratio. Set to 0 to skip computing sigma and use provided value for sig.
    int npts_acc = 1e4;     // number of points to estimate acceptance ratio. Not used if acc_target =0.
    int nburn = 1e2;        // number of burn-in points.  Not used if acc_target =0.



    /* ------------------------------------------------------ */
    /*   Code you might not want to change                    */
    /*     (except to change output printed to screen)        */
    /* ------------------------------------------------------ */

    // Projection string
    string projstr;
    if(projmethod == 1)  projstr = "Newton";
    else projstr = "Symmetric";

    cout << " ------------- " << endl;
    cout << "Beginning sampling!" << endl;
    cout << "example = " << example << endl;
    cout << "n = " << n << endl;    
    cout << "projection method = " << projstr << endl;


    /* Set up Equations object */
    Equations myeqns;  // this class holds all information needed by sampler. It will be populated by information in T. 
    T poly;
    VectorXd x0;
    x0.resize(n*poly.d);
    if(example == 1 || example == 2 || example == 3) {
        x0 = poly.default_initial(n);   // initializes example: initial condition, all parameters in class
    }
    else {
        x0 = poly.default_initial(n,iseed);   // initializes example: initial condition, all parameters in class
    }


    // initalize equations
    myeqns.nvars = poly.n*poly.d;         // # of variables
    myeqns.neqns = poly.m;                // # of equations
    myeqns.eval_q = std::bind(&T::q, poly,_1, _2);   // function which evaluates equations
    myeqns.eval_Dq = std::bind(&T::dq, poly,_1,_2);  // function which evaluates Jacobian of equations

    // Optionally bind to other functions
    if(example == 2) {
        myeqns.energy = std::bind(&EnergyFunctions::lattice_corner,_1,5,(int) round(sqrt(n)));   // energy function
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
    // sets sigma within mysampler object
    if(acc_target > 1e-6) {
        cout << "\nFinding sigma for target acceptance a = " << acc_target << " ..."<< endl;
        double sig_high = 1;
        mysampler.sig_finder(acc_target, sig_high, npts_acc, nburn);
    }
    


    /* Sample! */
    cout << "\nSampling...." << endl;
    mysampler.sample(npts,dsave);

    


    /* Write diagnostics to screen */
    cout << "Done sampling! " << endl;
    cout << "sig = " << mysampler.sig << endl;
    cout << "  This run took " << mysampler.time << " seconds" << endl;
    cout << "  Average acceptance ratio = " << mysampler.rej.acc_avg() << endl;
    cout << "  Average no. iterations in Projection method = " << mysampler.niteravg() << endl;
    cout << "     Avg no. iterations (successful)    =  "  << mysampler.niteravg(1) << endl;
    cout << "     Avg no. iterations (unsuccessful)  =  "  << mysampler.niteravg(-1) << endl;
    cout << "  sparsities (# nonzeros, Dq, Dq'*Dq, L): " << mysampler.sparsity().transpose() << endl;
    

    // Write rejection statistics to screen
    mysampler.print_rej();

    // Writing timing information to screen
    mysampler.print_time();   

    // Write statistics to file
    if(datafile != "") {
        cout << "Saving data to file" << endl;
        mysampler.write_stats(datafile);
    }
    

    // How to print out / access internal data (uncomment only when npts is small!)
    //cout << "  stats = \n" << mysampler.stats << endl;  // statistics
    //cout << "  niter = " << mysampler.niter.transpose() << endl;  // newton iterations (<0 when failed to converge)

    // write edges to file, for polygon
    if(example == 4) {
        ofstream myfile1 ("edges.txt");
        for(int i=0; i<poly.edges.rows(); i++) {
         myfile1 << poly.edges.row(i) << endl;
        }
         myfile1.close();
    }

    // print out matrices generated
    /*MatrixXd stats = mysampler.stats;
    MatrixXd M(poly.sl,poly.sl);
    for(int i=0; i<stats.rows(); i++) {
        M = poly.mtx_of_vec(stats.row(i));
        cout << "inner prod " << i << ":\n" << M.transpose() * M << endl;
        //cout << "stat " << i << ":\n" << poly.mtx_of_vec(stats.row(i)) << endl;
    }
    */

}

//}