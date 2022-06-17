// 
// Example of how to use ManifoldSampler
//
// This example sets up a quadrilateral with a triangle attached, in 2d space
// The quarilateral has vertices 0,1,2,3, and the triangle has vertices 2,4,5. 
// Avoid choosing a square, as its configuration space has a singularity that makes it hard to sample
// 
// Created April 9, 2021
// modified June 16, 2022
// 
//
//

#include <iostream>
#include <iomanip>    // for setting precision
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Framework.hpp"
#include "ManifoldSampler.hpp"

using namespace Eigen;
using namespace std;
using namespace std::placeholders;  // for binding functions

#define PI  3.141592653589793238462643383279502884




// different ways to define your own energy function here, if desired
// already defined, in EnergyFunctions
double EnergyFunctions::energy(const VectorXd& x, const Framework& myframe) { 
    return 1.0; 
}

// write your own function, with the right type
double myenergy(const VectorXd& x, const Framework& myframe) { 
    return 1.0; 
}



// Main loop
int main(int argc,char *argv[])
{   

    /* Set up the framework */

    // Parameters
    int d = 2;      // dimension of ambient space
    int n = 6;      // number of vertices
    int m = 7;      // number of edges
    MatrixXi edges;  // m x 2 matrix of edges
    VectorXd x0;     // initial condition, a dn-dimensional vector

    // Construct the edges
    edges.resize(m,2);
    edges << 0,1,  // square's edges
             1,2,
             2,3,
             3,0,
             2,4,  // triangle's edges
             2,5,
             4,5;


    // Construct the initial condition
    x0.resize(n*d);
    double th = PI/3;
    x0 << 0,0,  // square's vertices
          1,0,
          1+cos(th),sin(th),
          0,1,
          2,1,  // triangle's additional vertices
          1.5,1.5;

    // Create a framework object
    Framework myframework(n,d,edges);
    myframework.set_lengths_from_x(x0);   // set the lengths from given initial condition
    // If you want to set the lengths by hand, do this as:
    //      Framework myframework(n,d,edges,lengths);
    // where lengths is a m-dimensional VectorXd

    // Let's print out what the lengths are
    cout << "You constructed a Framework object with lengths =\n  " << myframework.lengths.transpose() << endl;


    /* Set up the sampler */

    // Sampler parameters
    int npts = 5e2;       // how many points to generate
    double sig = 0.3;    // step size
    double tol = 1e-5;    // tolerance for Newton's method   (optional)
    int maxIter = 40;     // max no. of iterations for Newton's method (optional)
    int dsave = 1;        // how often to save statistics (optional; default = 1)


    // Create a sampler object
    ManifoldSampler mysampler(myframework,x0,sig,tol,maxIter);


    // print sparsities
    cout << "sparsities: " << mysampler.sparsity().transpose() << endl;

    // Change the statistics that are saved
    mysampler.statsfcn = &StatisticsFunctions::pts;  // saves all points


    // Some optional things to change
    // Change projection method
          mysampler.projmethod = ManifoldSampler::cProjNewton;
          //mysampler.projmethod = ManifoldSampler::cProjSym;

    // Set the seed to a particular value, for reproducibility. Default device noise.
          //mysampler.setSeed(10012);   // here the seed is set to 10012.

    // include debug comments, or not (default = No)
          //mysampler.ifdebug = ManifoldSampler::cYes; 

    // Change energy function, if desired. 
    // Two examples here; one with no parameters, another with parameters, that must be bound to have the right type.
          mysampler.ffcn = &myenergy;  // a function defined in code above 
          //mysampler.ffcn = &EnergyFunctions::energy;  // a function defined in code above (could be below)
          //double kspring = 10;
          //mysampler.ffcn = bind(EnergyFunctions::bendingPolymer, _1, _2, kspring); // bending energy. You can write your own energy function. 

   
    // Sample!
    mysampler.sample(npts,dsave);


    // Write diagnostics to screen
    cout << "Done! " << endl;
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
    mysampler.write_stats("example1.txt");

    // How to print out / access internal data (uncomment only when npts is small!)
    //cout << "  stats = \n" << mysampler.stats << endl;  // statistics
    cout << "  niter = " << mysampler.niter.transpose() << endl;  // newton iterations (<0 when failed to converge)


}





