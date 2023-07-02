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



using namespace Eigen;
using namespace std;
using namespace std::placeholders;  // for binding functions

#define PI  3.141592653589793238462643383279502884


class Ellipsoid {
public:
    int neqns = 1;
    int dim = 3;      // dimension of ellipsoid (if you change, you should also change the axes vector)
    vector<double> a {1,1,3};      // axes of ellipsoid
    void q(const VectorXd& x, VectorXd& qvals) {    // constraint
        qvals.resize(1);
        qvals(0) = -1;
        for(int i=0; i<dim; i++) {
            qvals(0) += x(i)*x(i) / (a[i]*a[i]);
        }
    }
    void dq(const VectorXd& x, SpMat& Dq) {    // jacobian of constraint
        std::vector<Trip> tripletList;  // holds row,col,val for constructing jacobian
        tripletList.reserve(dim);
        for(int i=0; i<dim; i++) {
            tripletList.push_back(Trip(i,0, 2*x(i)/(a[i]*a[i])));
        }
        Dq.resize(dim,neqns);
        Dq.setFromTriplets(tripletList.begin(), tripletList.end());
    }
};



// Main loop
int main(int argc,char *argv[])
{   

    /* Set up Equations */
    Ellipsoid ellipse;
    Equations myeqns;
    myeqns.nvars = ellipse.dim;
    myeqns.neqns = ellipse.neqns;
    myeqns.eval_q = std::bind(&Ellipsoid::q, ellipse,_1, _2);
    myeqns.eval_Dq = std::bind(&Ellipsoid::dq, ellipse,_1,_2);

    VectorXd x0(3);
    x0 << ellipse.a[0],0,0;

 
    // Optionally bind to other functions
    //myeqns.energy = std::bind(&EnergyFunctions::unit_fcn,_1);
    //myeqns.evalstats = std::bind(&StatisticsFunctions::nostats,_1);  // turn off statistics; don't save anything


    /* Set up Sampler */
    // Sampler parameters
    int npts = 1e3;       // how many points to generate
    double sig = 0.7;    // step size
    double tol = 1e-5;    // tolerance for Newton's method   (optional)
    int maxIter = 100;     // max no. of iterations for Newton's method (optional)
    int dsave = 2;        // how often to save statistics (optional; default = 1)


    /* Create a sampler object */
    ManifoldSampler mysampler(myeqns,x0,sig,tol,maxIter);


    /* Some optional things to change */
    // Change projection method
    //mysampler.projmethod = ManifoldSampler::cProjNewton;
    //mysampler.projmethod = ManifoldSampler::cProjSym;    // this one is better

    // Set the seed to a particular value, for reproducibility. Default device noise.
    //mysampler.setSeed(10012);   // here the seed is set to 10012.

    // include debug comments, or not (default = No)
    //mysampler.ifdebug = ManifoldSampler::cYes; 

    // turn off the delta-function measure |grad q|^{-1}
    mysampler.ifqdet = ManifoldSampler::cNo;
   

    /* Sample! */
    mysampler.sample(npts,dsave);


    /* Write diagnostics to screen */
    cout << "Done! " << endl;
    cout << "sparsities: " << mysampler.sparsity().transpose() << endl;
    cout << "  This run took " << mysampler.time << " seconds" << endl;
    cout << "  Average acceptance ratio = " << mysampler.rej.acc_avg() << endl;
    cout << "  Average no. iterations in Projection method = " 
         << mysampler.niteravg() << endl;
    cout << "     Avg no. iterations (successful)    =  "  << mysampler.niteravg(1) << endl;
    cout << "     Avg no. iterations (unsuccessful)  =  "  << mysampler.niteravg(-1) << endl;
    

    // Write rejection statistics to screen
    mysampler.print_rej();

    // Writing timing information to screen
    mysampler.print_time();   

    // Write statistics to file
    mysampler.write_stats("Data/example2.txt");

    // How to print out / access internal data (uncomment only when npts is small!)
    //cout << "  stats = \n" << mysampler.stats << endl;  // statistics
    cout << "  niter = " << mysampler.niter.transpose() << endl;  // newton iterations (<0 when failed to converge)


}





