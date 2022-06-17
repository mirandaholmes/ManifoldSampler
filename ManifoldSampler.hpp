//
//  ManifoldSampler.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on June 16, 2022.
//
// 
//
//
//

#ifndef ManifoldSampler_hpp
#define ManifoldSampler_hpp

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>
#include <ctime>             // for timing operations
#include <random>            // generate random variables 
#include <iomanip>           // for setw
#include <fstream>           // for writing to files
#include <chrono>   // for high-precision timing
#include "Framework.hpp"

using namespace std;
using namespace Eigen;

typedef std::chrono::high_resolution_clock Clock;

// column ordering, for factoring sparse matrices
typedef AMDOrdering<int> MyOrdering; 
typedef COLAMDOrdering<int> COLAMD;   // this one is the cheapest (but doesn't work with Cholesky decomposition)

// Solvers for Newton step. 
typedef SparseLU<SpMat, COLAMD> NewtonFac;    // Solver for Newton step (faster for large n)

// Decompositions for linear equations in Newton solve. 
typedef SimplicialLLT<SpMat, Eigen::Lower, MyOrdering> Cholesky;  // Cholesky decomposition; for computing tangent steps
//typedef SimplicialLDLT<SpMat, Eigen::Lower, MyOrdering> Cholesky;  // Cholesky decomposition; alternative version


// aliases for energy functions and statistics function
using EnergyFcn = std::function<double(const VectorXd&, const Framework&)>; // type alias to std::function
using StatsFcn  = std::function<VectorXd(const VectorXd&, const Framework&)>; // type alias to std::function


// Class to hold rejection statistics for different kinds of rejections
class Rejections {
public:
	// Different kinds of rejections
    int project {0};
    int ineq {0};
    int metropolis {0};
    int revproj {0};
    int revval {0};
    int acc {0};

    // Total number of rejections
    int total(void) { return project + metropolis + revproj + revval + ineq;};
    // Acceptance rate
    double acc_avg(void) { return (double)acc / (double)npts(); };
    // Total number of points
    int npts(void) { return project + metropolis + revproj + revval + ineq + acc;};
    // Reset rejections
    void reset(void) { project=0; metropolis=0; revproj=0; revval=0; ineq=0;};

    // Write to screen
    void print(void);

    // Constructor
    Rejections(void) {} // default constructor
};



// Struct to hold Energy functions, so the user doesn't have to rewrite them
// All of these should have type "EnergyFcn", 
// OR have extra parameters so that using "bind" we obtain an EnergyFcn. 
// Each energy function returns an "energy" U; the probability we wish to sample is exp(-U). 
struct EnergyFunctions {
      static double unit_fcn(const VectorXd& x, const Framework& myframe) { return 0.0; }
      static double bendingPolymer(const VectorXd& x, const Framework& myframe, const double springconst);
      static double energy(const VectorXd& x, const Framework& myframe);   // use this to define your own energy function
};

// Struct to hold various statistics functions, so the user doesn't have to rewrite them
// All of these should have type "StatsFcn", 
// OR have extra parameters so that using "bind" we obtain a StatsFcn.
struct StatisticsFunctions {
    static VectorXd pts(const VectorXd& x, const Framework& myframe) { return x; };  // default statistics; =points
    static VectorXd nostats(const VectorXd& x, const Framework& frame) { 
            VectorXd data; data.resize(0); return data; }
};


class ManifoldSampler {
public:

    // Sampling parameters
    Framework frame;    // framework to sample
    VectorXd x0;        // initial point
    double sig;         // step size in tangent direction
    double tol;         // tolerance for constraints holding
    int maxIter;        // max # of iterations in newton's method
    double errfac = 0.95;      // how much qerr must decrease by each iteration of Newton / Symmetric Newton
    int seed = 0;       // seed for rng; 0 = based on device noise
    int projmethod = cProjSym;  // which projection method
    EnergyFcn ffcn { &EnergyFunctions::unit_fcn } ;  // energy function to sample on framework (default = 0)
    StatsFcn statsfcn { &StatisticsFunctions::pts }; // statistics to compute on framework (default = x)


    // Constructor
    ManifoldSampler(Framework frame0, 
                   VectorXd x0,
    	           double sigtan0, 
                   double tol0 = 1e-6,
 	               int maxIter0 = 40);

    // Sampling functions
    int sample(int npts=1, int dsave=1);   // the main sampling function
    int project_newton(const VectorXd&, const SpMat&, VectorXd&, SpMat&, NewtonFac&);  // project to manifold
    int project_symmetric(const VectorXd&, const SpMat&, VectorXd&, Cholesky&);  // doesn't update Jacobian
    
    // Functions to find an initial point on the manifold
    int find_point_on_manifold(const VectorXd&, VectorXd&);  // tries to output a point on the manifold
    void steepestDescent(const VectorXd&, VectorXd&, double s=1e-4, int nsteps = 1);  // steepest descent for a given number of steps 


    // Things to keep track of during sampling
    Rejections rej;   // class holding rejection statistics
    MatrixXd stats;   // statistics, set by a user-defined function
    int nstats {0};   // how many statistics to keep track of 
    int nsave;        // total number of saved points
    double time;      // total time the for loop took
    VectorXi niter;   // number of Newton iterations at each step
    int ifSaveNIter {cYes};  // save list of niter values

    // Other options for running the code
    void  setSeed(int);      // set the rng seed
    int ifdebug   {cNo};
    int prec_pts   = 6;      // precision for writing output
    int write_stats(string statsfile, int dwrite=1);  // write statistics to file
    void print_time(void);        // print timing info to screen
    void print_rej(void);         // print rejections info to screen
    VectorXd sparsity(void);    // return sparsities of Dq, Dq'*Dq, L
    double niteravg(int opt=0);   // print average niter (option to print when successful, or unsuccessful)

    
    // Generate random numbers
    void   initializeRandom(void);   // initialize random number generators (use after change seed)
    double unifRand(void);           // output a uniform random number on [0,1]
    double normRand(void);           // output a standard normal random number


    // Helper functions
    void cholesky_diagonal(const SimplicialLLT<SpMat, Eigen::Lower, MyOrdering>&, VectorXd&); // extracts diagonal of cholesky matrix
    void cholesky_diagonal(const SimplicialLDLT<SpMat, Eigen::Lower, MyOrdering>&, VectorXd&); // extracts diagonal of cholesky matrix

    // Flags
    static const int cSuccess { 0 };  // test delete -- maybe not needed
    static const int cFail    { -1 };
    static const int cYes     { 0 };
    static const int cNo    { -1 };
    static const int cProjNewton { 0 };
    static const int cProjSym  { 1 };


private:     
    // For generating random numbers
    mt19937                       rng;         // defines the generator
    normal_distribution<double>   mZdist;      // standard random numbers
    uniform_real_distribution<double> mUdist;  // uniform on [0,1]
};


#endif /* ManifoldSampler_hpp */

