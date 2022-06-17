//
//  Framework.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on June 16, 2022.
//
// 
//
//
//


#ifndef Framework_hpp
#define Framework_hpp

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>


using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;       // for constructing sparse matrices


class Framework {
public:
	    
    const int n;   // number of vertices
    const int d;   // dimension of vertices
    const int m;   // number of edges/constraints
    const MatrixXi pairs;  // mx2 list of pairs of vertices that are connected by edges
    VectorXd lengths;  // lengths of edges


    void eval_q (const VectorXd&, VectorXd&);  // evaluate constraints 
    void eval_Dq(const VectorXd&, SpMat&);     // evaluate Jacobian of constraints 

    void set_lengths_from_x(const VectorXd& x);  // set lengths to those from current configuration x

    // Constructors
    Framework(int, int, MatrixXi&, VectorXd&);  // Set n, d, pairs, lengths (in that order)
    Framework(int, int, MatrixXi&);             // Set n, d, pairs
};


#endif /* Framework_hpp */


