//
//  Framework.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on June 16, 2022.
//
// 
//
//
//

#include "Framework.hpp"

// Constructor #1: with lengths
Framework::Framework(int n0, int d0, MatrixXi& pairs0, VectorXd& lengths0)
: n(n0), d(d0), pairs(pairs0), lengths(lengths0), m(pairs0.rows()) {
}

// Constructor #2: no lengths; initialize lengths to 0
Framework::Framework(int n0, int d0, MatrixXi& pairs0)
: n(n0), d(d0), pairs(pairs0), m(pairs0.rows()) {
    lengths = VectorXd::Zero(m);
}



// Set lengths to the ones obtained from current configuration x
void Framework::set_lengths_from_x(const VectorXd& x) {
    VectorXd q(m); 
    eval_q(x,q);
    lengths = q.cwiseSqrt();
}


// Evaluate constraints
void Framework::eval_q(const VectorXd& x, VectorXd& q) {
	int i,j;
	double a, d2;
    q.resize(m);
    for (int k=0; k < m; k++) {
        i = pairs(k,0);
        j = pairs(k,1);
        a = lengths(k);
        q(k) = -a*a;
        for (int l=0; l<d; l++) {
        	d2 = x(i*d+l) - x(j*d+l);
        	q(k) += d2*d2;
        }
    }
}


// Evaluate Jacobian of constraints, Dq = \grad q. 
// Columns of Dq are gradients of each constraints. 
// Jacobian matrix is sparse, constructed using a vector of triplets, each constructed as 
//   Trip(row,col,val). 
//   A detailed example is here: http://eigen.tuxfamily.org/dox/TutorialSparse_example_details.html
//
// Notes: 
// * It might be faster to pre-allocate nonzeros per column, then fill in element-by-element, 
//   as in example 2 on this page: http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling
//   Could do this in the constructor. 
// * Alternatively, we could construct Dq initially, extract index order, then iterate through 
//   that natural index order to evaluate & set values, as in this example: http://eigen.tuxfamily.org/dox/group__TutorialSparse.html
//  (see "iterating over nonzero coefficients")
//   
void Framework::eval_Dq(const VectorXd& x, SpMat& Dq) {
	std::vector<Trip> tripletList;  // holds row,col,val for constructing jacobian
	tripletList.reserve(m*(2*d));
	int i,j;
	double val;
	for (int k=0; k < m; k++) {
        i = pairs(k,0);
        j = pairs(k,1);
        for (int l=0; l<d; l++) {
            val = x(i*d+l) - x(j*d+l);
            tripletList.push_back(Trip(i*d+l,k,2*val));
            tripletList.push_back(Trip(j*d+l,k,-2*val));
        }
    }
    Dq.resize(d*n,m);
    Dq.setFromTriplets(tripletList.begin(), tripletList.end());
}
