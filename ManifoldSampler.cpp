//
//  ManifoldSampler.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on June 16, 2022.
//
// 
//
//
//

#include "ManifoldSampler.hpp"



// Constructor
ManifoldSampler::ManifoldSampler(Equations& eqns0, 
                   VectorXd x0,
    	           double sig0, 
                   double tol0,
 	               int maxIter0)
: eqns(eqns0), x0(x0), sig(sig0), 
  tol(tol0), maxIter(maxIter0), 
  mZdist(0.,1.), mUdist(0.0,1.0)
{
    initializeRandom();     // initialize random number generators
}


// ====================================================
//           Main Sampling Algorithm
// ====================================================
//
// input: 
//      npts = number of points to sample (default = 1)
//      dsave = how often to save statistics (default = 1)
// 
// returns: 
//      avg acceptance rate
//
// saves:
//      x0 = last point sampled
//
// 
double ManifoldSampler::sample(int npts, int dsave) {

	// extract equations info, for ease of reference
	int m = eqns.neqns;    // number of constraints
	int nvars = eqns.nvars;    // number of variables

    // Variables needed in code
    SpMat Dq1, Dq2, Dqxrev;              // Jacobians, at x, y, and for reverse move
    SpMat* pDqx = &Dq1;                  // pointer to jacobian, x
    SpMat* pDqy = &Dq2;                  // pointer to jacobian, y
    SpMat* pDqtmp;                       // temporary; used for switching
	VectorXd x(nvars), y(nvars), xrev(nvars);   // points on manifold
    VectorXd r(nvars);                   // full random step
    VectorXd z1,z2;                      // help in solving for tangent steps
    VectorXd vtan(nvars), vtany(nvars);  // tangent steps, forward & reverse
    double acc;                          // total metropolis ratio
    double vdiff,ux,uy,udiff, qdet_yx;   // factors used to compute metropolis ratio
    VectorXd diagx(m), diagy(m);         // diagonals of cholesky matrix (for computing det(Dqy)/det(Dqx))
    int flag;                            // success flag
    clock_t start_time, end_time;        // start,end times


    // Construct initial Dqx, energy
    x = x0;
    eqns.eval_Dq(x,Dq1);
    ux = eqns.energy(x);


    // Construct initial Cholesky objects, and perform initial analysis 
    Cholesky chol1, chol2;   // cholesky desompositions of Dq^T*Dq, in x and y
    Cholesky* pcholx = &chol1;  // pointer to cholesky object, x
    Cholesky* pcholy = &chol2;  // pointer to cholesky object, y
    Cholesky* pcholtmp;         // temporary; used for switching
    chol1.analyzePattern(Dq1.transpose()*Dq1);    // depends only on sparsity pattern, not numerical values
    chol2.analyzePattern(Dq1.transpose()*Dq1);

    // Initial Cholesky factorization in x
    (*pcholx).factorize(Dq1.transpose()*Dq1);    // factorize these using numerical values, to start
    cholesky_diagonal((*pcholx),diagx);   // sets diagx. Replaces: qx = sqrt(chol1.determinant());  


    // Object used in Newton's method to solve linear equations
    // set here so its sparsity pattern is analyzed beforehand
    NewtonFac  newtonfac;    
    if(projmethod == cProjNewton) {
        newtonfac.analyzePattern(Dq1.transpose()*Dq1);
    }  
    

    // Initialize matrices to hold data
    nsave = floor(npts / dsave) + 1;
    nstats = eqns.evalstats(x).size();   // how many statistics to calculate
    if(npts%dsave==0) nsave--; 
    stats = MatrixXd::Zero(nsave,nstats);
    int ipt = 0;
    if(ifSaveNIter == cYes) niter.resize(npts,1);
    rej.reset();   // reset rejection statistics (could remove this line, if you want to accumulate from run to run)


    // *************  Loop through MCMC steps  ***************
    start_time = clock();            // start timing now
    for(int i=0; i<npts; i++) {

    	// =============  Save data  =============
        if(i%dsave==0 && nstats >0) {
            stats.row(ipt) = eqns.evalstats(x);
            ipt++;
        }

    	// =============  Get tangent steps  =============
    	for (int j=0; j<nvars; j++) r(j) = normRand(); // total random step
        z1 = (*pcholx).solve((*pDqx).transpose()*r);
        vtan = r - (*pDqx)*z1;
        vtan = sig*vtan;     // scale to have correct magnitude



        // TEST -- compare ninter for different projection methods (uncomment to test)
        //cout << setw(6) << project_newton(x+vtan,*pDqx,y,*pDqy,newtonfac);
        //cout << setw(6) << project_symmetric(x+vtan,*pDqx,y,*pcholx) << endl;


    	// =============  Take step & project  =============
        if(projmethod == cProjNewton) {
            flag = project_newton(x+vtan,*pDqx,y,*pDqy,newtonfac);  // sets y,Dqy
        }
        if(projmethod == cProjSym) {
            flag = project_symmetric(x+vtan,*pDqx,y,*pcholx);
        }

    	if(ifSaveNIter == cYes) niter(i) = flag;  // save number of iterations. negative if failed to converge.

    	if(flag < 0) {  // projection failed to converge to a solution
    	    rej.project++;
    	    continue;
        }
    	
    	// DEBUG
        if(ifdebug==cYes) {
        	cout << "iter:" << i << endl;
    		//cout << " y - z = " << (y - x - vtan).norm() << endl;
    	}


    	// =============  Inequality Check  =============
    	// **** (no inequalities for now) ******
    	flag = cSuccess;
    	if(flag == cFail) {
    		rej.ineq++;
    		continue;
    	}


        // =============  Get reverse tangent step  =============
        if(projmethod == cProjSym) eqns.eval_Dq(y,*pDqy);  // evaluate Dqy
        (*pcholy).factorize((*pDqy).transpose()*(*pDqy));    // cholesky decomposition, with values
        r = x-y;
        z1 = (*pcholy).solve((*pDqy).transpose()*r);
        vtany = r - (*pDqy)*z1;


    	// =============  Metropolis step  =============
        vdiff = 0.5*( vtany.squaredNorm()/(sig*sig) - vtan.squaredNorm()/(sig*sig));
        uy = eqns.energy(y);
        udiff = uy-ux;   // energy difference
        cholesky_diagonal((*pcholy), diagy);  // sets diagy. Replaces: qy = sqrt(chol2.determinant());  
    	qdet_yx = (diagy.cwiseQuotient(diagx)).array().prod();
        acc = exp(-vdiff) * exp(-udiff) * 1./(qdet_yx);  // metropolis acceptance ratio


        // DEBUG
    	if(ifdebug == cYes) {
    		cout << "  acc = " << acc << endl;
            cout << "  uy = " << uy << endl;
            cout << "  ux = " << ux << endl;
    		cout << "  fy/fx=exp(-udiff) = " << exp(-(uy-ux)) << endl;
    		cout << "  vy/vx       = " << exp(-vdiff) << endl;
            cout << "    prev      = " << exp(-vtany.squaredNorm()/(2*sig*sig))/exp(-vtan.squaredNorm()/(2*sig*sig)) << endl;
            cout << "  qdet_yx     = " << qdet_yx << endl;
    		cout << "  |vtan|^2/|vtany|^2   = " << vtan.squaredNorm()/vtany.squaredNorm() << endl;
        }
		

        if(unifRand() > acc) {
        	rej.metropolis++;
            continue;  
        }


    	// =============  Reverse projection  =============
        if(projmethod == cProjNewton) {
            flag = project_newton(y+vtany,*pDqy,xrev,Dqxrev,newtonfac);  
        }
    	if(projmethod == cProjSym) {
            flag = project_symmetric(y+vtany,*pDqy,xrev,*pcholy); 
        }

    	if(flag < 0) {  // projection failed to converge
    		rej.revproj++;
    		continue;
    	}
    	if((x-xrev).norm() > tol*nvars*10) {  // projection converges, but to incorrect x (10 is arbitrary)
            rej.revval++;
            continue;
    	}
        

    	// =============  Accept move; copy over data to x  ============= 
        rej.acc++;
    	x = y;
        // switch jacobian pointers
        pDqtmp = pDqx;    
        pDqx = pDqy;
        pDqy = pDqtmp;    
    	ux = uy;
    	diagx = diagy;   // Replaces: qx = qy;
        // switch cholesky pointers
    	pcholtmp = pcholx;
        pcholx = pcholy;
        pcholy = pcholtmp;
    }

    // save time
    time = (clock() - start_time ) / (double) CLOCKS_PER_SEC;

    // save last point
    x0 = x;

    // return avg acceptance rate
    return rej.acc_avg();
}


// ====================================================
//           Project to manifold
// ====================================================
//
// Solve for variables a such that q(z + Qx*a) = 0. 
//

// project_newton: 
//    Runs Traditional Newton's method, to solve system of equations
// Sets:
//        y = z+Qx*a, 
//        Dqy 
// Returns: 
//        iter if successful; -iter if unsuccessful, where iter = # of iterations
// NewtonFac must have run analyzePattern() already. 
// 
int ManifoldSampler::project_newton(const VectorXd& z, const SpMat& Qx, VectorXd& y, SpMat& Qy, NewtonFac& projsolver) {

    int nvars = Qx.rows();
    int neqns = Qx.cols();
    VectorXd q(neqns), da(neqns), a(neqns);
    double qerr, qerr0;

    // Initial guess; all 0's
    a.setZero();   
    y = z + Qx*a;
    eqns.eval_Dq(y,Qy);
    qerr = std::numeric_limits<double>::infinity();;

    // Loop until get close enough to a solution
    for(int iter=0; iter < maxIter; iter++) {
        // evaluate q
        eqns.eval_q(y,q);
        qerr0 = qerr;       // save old value, to test if it decreases
        qerr = q.cwiseAbs().maxCoeff();  // maximum error, pointwise

        // Check if converged
        if(qerr < tol) {
            return iter;  // y, Qy are already up to date
        }

        // Check if qerr failed to decrease sufficiently
        if(qerr > errfac*qerr0) {
            return -iter;
        }

        // Otherwise, evaluate Jacobian and solve for increment
        eqns.eval_Dq(y,Qy);
        projsolver.factorize(Qy.transpose()*Qx);
        da = projsolver.solve(-q);
        a = a + da;
        y = z + Qx*a;
    }

    return -maxIter;  // failed to find a solution 
}


// project_symmetric: 
//    Runs Symmetric Newton's method, to solve system of equations
// Sets:
//        y = z+Qx*a, 
// Returns: 
//        iter if successful; -iter if unsuccessful, where iter = # of iterations
// Cholesky must have run analyzePattern() already. 
// 
int ManifoldSampler::project_symmetric(const VectorXd& z, const SpMat& Qx, VectorXd& y, Cholesky& projsolver) {

    int nvars = Qx.rows();
    int neqns = Qx.cols();
    VectorXd q(neqns), da(neqns), a(neqns);
    double qerr, qerr0;

    // Initial guess; all 0's
    a.setZero();   
    y = z + Qx*a;
    qerr = std::numeric_limits<double>::infinity();;

    // Loop until get close enough to a solution
    for(int iter=0; iter < maxIter; iter++) {
        // evaluate q
        eqns.eval_q(y,q);
        qerr0 = qerr;       // save old value, to test if it decreases
        qerr = q.cwiseAbs().maxCoeff();  // maximum error, pointwise

        // Check if converged
        if(qerr < tol) {
            return iter;  // y is already up to date
        }

        // Check if qerr failed to decrease sufficiently
        if(qerr > errfac*qerr0) {
            return -iter;
        }

        // Otherwise, solve for increment
        da = projsolver.solve(-q);
        a = a + da;
        y = z + Qx*a;
    }

    return -maxIter;  // failed to find a solution
}




// ====================================================
//         Construct initial point on manifold
// ====================================================

// Try to find a point x on the manifold, using Newton's method
// Solves system of equations for vector a:
//    q(x0 + Dq(x0)*a) = 0
// Output is x = x0 + Dq(x0)*a.
// Returns flag indicating if projection succeeded.
// 
int ManifoldSampler::find_point_on_manifold(const VectorXd& x0, VectorXd& x) {
    int flag;
    SpMat Dq,Dqy;
    eqns.eval_Dq(x0,Dq);
    NewtonFac projsolver;
    projsolver.analyzePattern(Dq.transpose()*Dq);
    flag = project_newton(x0, Dq, x, Dqy, projsolver);
    return flag;
}

// Do a few steps of steepest descent, minimizing the energy 
//    E = 0.5*\sum_i q_i^2. 
// Use this if find_point_on_manifold fails. 
//
void ManifoldSampler::steepestDescent(const VectorXd& x0, VectorXd& x, double s, int nsteps) {
    SpMat Dq;
    VectorXd q(eqns.neqns);
    x = x0;
    for(int i=0; i<nsteps; i++) {
        eqns.eval_q(x,q);
        eqns.eval_Dq(x,Dq);
        x = x - s*Dq*q;
    }
}


// ====================================================
//           Helper Linear Algebra functions
// ====================================================

// extracts diagonal of Cholesky L-matrix
void ManifoldSampler::cholesky_diagonal(const SimplicialLLT<SpMat, Eigen::Lower, MyOrdering>& chol, VectorXd& diag) {
    SparseMatrix<double,RowMajor> Lr = chol.matrixL();  // row-major, sorted
    SpMat L = Lr;  // column-major
    int ct = 0;
    for (int k=0; k<L.outerSize(); ++k) {
        for (SpMat::InnerIterator it(L,k); it; ++it) {
            if(it.row() == it.col()) {
                diag(ct) = it.value();
                ct++;
            }
        }
    }
}

// extracts diagonal of Cholesky L-matrix
void ManifoldSampler::cholesky_diagonal(const SimplicialLDLT<SpMat, Eigen::Lower, MyOrdering>& chol, VectorXd& diag) {
    diag = chol.vectorD().cwiseSqrt();
}





// ====================================================
//           Output Functions
// ====================================================

VectorXd ManifoldSampler::sparsity(void) {
    SpMat Dq;
    Cholesky chol;
    eqns.eval_Dq(x0,Dq);
    chol.analyzePattern(Dq.transpose()*Dq);
    chol.factorize(Dq.transpose()*Dq);
    SparseMatrix<double,RowMajor> Lr = chol.matrixL();  // row-major, sorted
    SpMat Dq2 = Dq.transpose()*Dq;

    VectorXd output(3);

    output(0) = Dq.nonZeros();
    output(1) = Dq2.nonZeros();
    output(2) = Lr.nonZeros();

    return output;
}


// compute average niter 
// opt = 0: average over all values
// opt = 1: average over successful values
// opt = -1: average over unsuccessful values
double ManifoldSampler::niteravg(int opt) {
    int ct = 0;
    double avg = 0;
    for(int i=0; i<niter.size(); i++) {
        if(opt == 0) {
            avg += abs(niter(i));
            ct++;
        }
        if(opt == 1 && niter(i) >= 0) {
            avg += niter(i);
            ct++;
        }
        if(opt == -1 && niter(i) <= 0) {
            avg += -niter(i);
            ct++;
        }
    }
    return avg/(double)ct;
}


// Write statistics to file, at intervals of dwrite points
int ManifoldSampler::write_stats(string statsfile,
                              int dwrite) {
    ofstream myfile1 (statsfile); // SHOULD check for problems
    myfile1.precision(prec_pts);   // Set precision

    // Write points to file
    cout << "write_stats: writing stats to file " << statsfile<< " at intervals of " << dwrite << " points." << endl;
    auto t1 = Clock::now();
    int i=0;
    while(i<nsave) {
        myfile1 << stats.row(i) << endl;
        i += dwrite;
    }
    myfile1.close();

    auto t2 = Clock::now();
    cout << "Writing to file took "
    << chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9 << " seconds" << std::endl;
    
    // return
    return ManifoldSampler::cSuccess;
}


// Print timing info to screen
void ManifoldSampler::print_time(void) {
	// Return to default precision
    cout << "-----  Timing  -----" << endl;
    cout << "Total time    : " << time << " seconds" << endl;
    //cout << "Time per point: " << mtime/mnpts << " seconds" << endl;
    double nn1 = 1e6;
    int npts = rej.npts();
    cout << "You can generate " << nn1 << " points in " << time*nn1/(double)npts << " seconds";
    cout << " or " << time*nn1/(double)npts/60.0 << " minutes." << endl;
}

// Print rejection statistics, by calling function in Rejections class
void ManifoldSampler::print_rej(void) {
	rej.print();
}

// Write rejection statistics to screen
void Rejections::print(void) {
	int w  = 12;   // width of boxes in table
	int w2 = 16;   // wider boxes

	// Set precision
    streamsize prec_init = std::cout.precision();  // current precision
    cout.precision(3);

    cout << "\n---------------  Rejection Statistics:  ---------------" << endl;
    cout << "Total number of points : " << npts()  << endl;
    cout << setw(w2) << "" << setw(w) << "Number" 
         << setw(w) << "Rate(Rej)" << setw(w) << "Rate(Tot)" << endl;
    cout << setw(w2) << "Project" << setw(w) << project 
         << setw(w) << (double)project/total() << setw(w) << (double)project/npts() << endl;
    cout << setw(w2) << "Ineq" << setw(w) << ineq
         << setw(w) << (double)ineq/total() << setw(w) << (double)ineq/npts() << endl;
    cout << setw(w2) << "Metropolis" << setw(w) << metropolis
         << setw(w) << (double)metropolis/total() << setw(w) << (double)metropolis/npts() << endl;
    cout << setw(w2) << "Reverse Project" << setw(w) << revproj
         << setw(w) << (double)revproj/total() << setw(w) << (double)revproj/npts() << endl;
    cout << setw(w2) << "Reverse Val" << setw(w) << revval
         << setw(w) << (double)revval/total() << setw(w) << (double)revval/npts() << endl;
    cout << setw(w2) << " --------" << endl;
    //cout << " " << endl;
    cout << setw(w2) << "Accept" << setw(w) << acc
         << setw(w) << ""<< setw(w) << (double)acc/npts() << endl;
    cout << " -------------------------------------------------------" << endl;

    // Return to default precision
    cout.precision(prec_init);
}


// ====================================================
//           Random number generators
// ====================================================


// Initialize random number generator
void  ManifoldSampler::initializeRandom(void) {
    if(seed==0) {
        seed = random_device()(); // choose seed randomly based on device configurations
    }
    rng.seed(seed);    // initialize generator for random number stream
}

// Output a uniform random number on [0,1]
double ManifoldSampler::unifRand(void) {
    return mUdist(rng);
}

// Output a standard normal random number 
double ManifoldSampler::normRand(void) {
    return mZdist(rng);
}

void ManifoldSampler::setSeed(int seed0) {
	seed = seed0;
	initializeRandom();
}





