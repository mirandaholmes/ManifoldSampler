# ManifoldSampler
This code implements a Markov Chain Monte Carlo algorithm to sample a probability distribution on a manifold
in a numerically efficient way. The manifold is defined by the level set of constraint functions: 
$$q(x)=0,$$
where $q:\mathbb{R}^n\to \mathbb{R}^m$, with $m < n$. The probability distribution to be sampled has the form
$$\pi(x) \propto f(x) |\nabla q|^{-1}\mu $$
(if the constraints are imposed by $\delta$-functions), where $\mu$ is the natural Hausdorff measure on the manifold, or 
$$\pi(x) \propto f(x) \mu $$
if the constraints are so-called "hard" constraints, that are imposed directly. (For the latter, be sure to set `ifqdet=cNo` in the ManifoldSampler class.)

The algorithm is described in [insert ref]. 

Two examples to show how to use the code are 
* example1.cpp: samples points uniformly on the surface of a d-dimensional ellipsoid (current version has d=3).
* example2.cpp: samples a simple framework: a collection of vertices with fixed-length edges.
Code to plot the output is in the Data/ folder. 

To compile this code, you will need to first install the Eigen library and add it to your compiler path. 
Eigen can be downloaded here: 
https://eigen.tuxfamily.org/index.php?title=Main_Page. 
It is recommended to compile with option -O3 for optimal efficiency. A makefile is included, but note that it includes more files than are strictly necessary, because it is set up to 
handle many examples. 

The examples described in the paper are found in codes in the PaperCodes/ folder. All examples can each be run using code
* examples_paper.cpp

To run a given example, make sure to set `typedef T` (lines 33-36) to the desired example, 
and set `int example=` to the corresponding integer (line 47). 


The data for the paper was generated by running
* run_polymer.cpp
* run_lattice.cpp
* run_matrix.cpp
* run_polygon.cpp

with appropriate parameters; these are also in the PaperCodes/ folder. 
The additional makefile in this folder can be used to compile all the paper-specific codes. 



The heart of the algorithm is in the ManifoldSampler class. 

ManifoldSampler must be constructed with an Equations object, an initial condition, and a step size parameter. 
All other parameters are optional and may be changed directly in main{}. 
The initial condition can be set in main{}, or, it can be 
constructed in a separate class file (as in option 3 below). 

The Equations class is a container to hold the equations defining the manifold, via member function `eval_q`, and the Jacobian 
of these equations ($\nabla q$), via member function `eval_Dq`. There are several options for how to fill in the content for this class. 
1. Write them directly in the Equations class. For example, you could create an Equations.cpp file
containing definitions of the functions and Jacobian. (There is no such example provided, because we wish 
to construct several examples in the same folder.)
2. Write them directly above main{}. (This is not very convenient as it is hard to vary parameters.)
3. Create another class containing function definitions, and bind the Equation class's functions
to this class's functions in main{}. All examples included do this: example2.cpp creates a new class directly above main{}, 
which other examples (example1.cpp, examples_paper.cpp, and all the run_[example].cpp files) have a separate class file, 
for readability. 

By default, the sampler samples a density with respect to the Hausdorff measure on the manifold which has  the form 
    $$\pi(x) \propto f(x) |\nabla q|^{-1},$$ 
where $f(x)$ is a function of type EnergyFcn, and $|\nabla q|^{-1}$ is the pseudodeterminant of the Jacobian of the constraints. 
This is the natural measure that arises when one considers physically-based sampling problems where constraints are 
approximations for stiff forces. If your constraints are so-called "hard" constraints, that don't come from delta-functions, 
then you don't need this pseudodeterminant factor. In this case, set `ifqdet=cNo` in ManifoldSampler.hpp. 
