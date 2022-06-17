// 
// Develop ManifoldSampler class
// 
// Created Jan 20, 2021
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


// Energy functions
namespace Trimer {
    double anglespring(const VectorXd& x, double kspring) {
        double costh = (x(0)-x(2))*(x(4)-x(2)) + (x(1)-x(3))*(x(5)-x(3));
        return kspring*(1-costh);   // wants to be flat
    }
}


// Main loop
int main(int argc,char *argv[])
{

    double th = PI/3;  // angle for polymer bending

    int d=2;   // dimension of polymer
    int n,m;
    double sig;
    MatrixXi edges;
    VectorXd x0;

    int numexpt = 1;
	VectorXd nlist(numexpt);   // list of n-values 
    VectorXd siglist(numexpt);  // list of sigmas  
    nlist << 3;        // list of sizes to consider; numexpt values
    siglist << 0.5;    // corresponding  sigmas; numexpt values

    // Sampler parameters
    int npts = 1e3;   // how many points to generate
    double tol = 1e-6;    // tolerance for Newton's method   (optional)
    int maxIter = 50;     // max no. of iterations for Newton's method (optional)
    int dsave = 1;     // how often to save statistics (optional; default = 1)

    double pert = 0;  // size of perturbation (random)


    // Stuff to keep track of during sampling
    VectorXd times(numexpt);
    VectorXd accept(numexpt); 
    VectorXd niter(numexpt);


    for (int i=0; i < numexpt; i++) {

        n = nlist(i);      // number of particles
        sig = siglist(i);  // sigma
        m = n-1;       // number of constraints
        cout << " ------ n = " << n << ", sig = " << sig << " ------ " << endl;

        // set edges
        edges.resize(m,2); edges.setZero();
        edges.col(0).setLinSpaced(m,0,m-1);
        edges.col(1).setLinSpaced(m,1,m);

        // construct basis vectors
        VectorXd e1 = VectorXd::Zero(d);
        VectorXd e2 = VectorXd::Zero(d);
        e1(0) = 1;
        e2(0) = cos(th);
        e2(1) = sin(th);

        // construct x0
        x0.resize(n*d);
        x0.setZero();
        x0.segment(d,d) = x0.head(d) + e2;
        for(int j=2; j<n; j++) {
            if(j%2 == 0) x0.segment(j*d,d) = x0.segment((j-1)*d,d) + e1;
            else x0.segment((j)*d,d) = x0.segment((j-1)*d,d) + e2;
        }

        // perturb, if desired
        //srand(seed);
        for(int j=0; j<n*d; j++) {
            x0(j) += pert* rand()/RAND_MAX;
        }


        // Create a framework object
        Framework myframework(n,d,edges);
        myframework.set_lengths_from_x(x0);   // set the lengths from given initial condition


        // Set up Equations 
        Equations myeqns;
        myeqns.nvars = myframework.nvars;
        myeqns.neqns = myframework.m;
        myeqns.eval_q = std::bind(&Framework::eval_q, myframework,_1,_2);
        myeqns.eval_Dq = std::bind(&Framework::eval_Dq, myframework,_1,_2);

        // Set energy function 
        double kspring = 7;
        myeqns.energy = std::bind(&Trimer::anglespring,_1,kspring);


        // Create a sampler object
        ManifoldSampler mysampler(myeqns,x0,sig,tol,maxIter);


        // Change projection method
        //mysampler.projmethod = ManifoldSampler::cProjNewton;
        mysampler.projmethod = ManifoldSampler::cProjSym;

        // Change seed
        //mysampler.setSeed(10012);

        // include debug comments, or not
        //mysampler.ifdebug = ManifoldSampler::cYes; 


        // Sample!
        mysampler.sample(npts,dsave);

        // Save timing
        times(i) = mysampler.time;
        accept(i) = mysampler.rej.acc_avg();
        niter(i) = mysampler.niteravg();
        

        // IN CASE INTERRUPTED
        cout << setprecision(3);
        cout << "times =  " << times.transpose() << endl;
        cout << "accept = " << accept.transpose() << endl;
        cout << "niter = " << niter.transpose() << endl;
        cout << setprecision(6);

        
        // print stats to screen
        mysampler.print_rej();
        mysampler.print_time();   

        cout << "  Average acceptance ratio = " << mysampler.rej.acc_avg() << endl;
        cout << "  Average no. iterations in Projection method = "  << mysampler.niteravg() << endl;
        cout << "     Avg no. iterations (successful)    =  "  << mysampler.niteravg(1) << endl;
        cout << "     Avg no. iterations (unsuccessful)  =  "  << mysampler.niteravg(-1) << endl;

        // How to print out / access internal data (uncomment only when npts is small!)
        //cout << "  stats = \n" << mysampler.stats << endl;  // statistics
        //cout << "  niter = " << mysampler.niter.transpose() << endl;  // newton iterations (<0 when failed to converge)

        // Write statistics to file
        mysampler.write_stats("Data/polymer.txt");

    }

    cout << setprecision(3);
    cout << "SUMMARY:" << endl;
    cout << "times =  " << times.transpose() << endl;
    cout << "accept = " << accept.transpose() << endl;
    cout << "niter = " << niter.transpose() << endl;


}
