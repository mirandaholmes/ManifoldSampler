// 
// Container to hold information about equations, needed by ManifoldSampler
// 
// You can write these functions directly, below the class definition, or alternatively
// you can bind these functions in another script, as in example1.cpp.
//
// 
// Created June 17, 2022 by Miranda Holmes-Cerfon
//
//
//
//

#ifndef Equations_hpp
#define Equations_hpp

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <functional> 


using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;       // for constructing sparse matrices


typedef std::function<void (const VectorXd&, VectorXd&)> Constraints;
typedef std::function<void (const VectorXd&, SpMat&)> Jacobian;
typedef std::function<double(const VectorXd&)> EnergyFcn; 
typedef std::function<VectorXd(const VectorXd&)> StatsFcn; 



// Struct to hold Energy functions, so the user doesn't have to rewrite them
// Should have type "EnergyFcn", OR have extra parameters so using "bind" gives a StatsFcn.
// Each energy function returns an "energy" U; the probability we wish to sample is exp(-U). 
struct EnergyFunctions {
	static double unit_fcn(const VectorXd& x) { return 0.0; }
};



// Struct to hold various statistics functions, so the user doesn't have to rewrite them
// Should have type "StatsFcn", OR have extra parameters so using "bind" gives a StatsFcn.
struct StatisticsFunctions {
    static VectorXd pts(const VectorXd& x) { return x; };  // default statistics; =points
    static VectorXd nostats(const VectorXd& x) { VectorXd data; data.resize(0); return data; }  // no statistics
};



class Equations {
public:
    int nvars;            // number of variables
    int neqns;            // number of constraints
    Constraints eval_q;   // system of equations to evaluate
    Jacobian eval_Dq;     // Jacobian (transpose) of system of equations
    EnergyFcn energy { &EnergyFunctions::unit_fcn };     // energy of system
    StatsFcn evalstats { &StatisticsFunctions::pts };    // statistics to evaluate
};


#endif /* Equations_hpp */

