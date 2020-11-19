#ifndef HEUNC_HPP_
#define HEUNC_HPP_

/* Class to numerically calculate the confluent Heun function for a real argument 
   Converted from MATLAB code
   Jamie Bamber 2020
  
   for the author of MATLAB code see below
*/

/* Comments from original author of MATLAB script. That is available here: https://github.com/motygin/confluent_Heun_functions */

// confluent Heun function, the first local solution of the equation
// HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
// at z=0 such that HeunC(0)=1 and
// HeunC'(0)=-q/gamma when gamma is not equal to 0
// HeunC'(z)/log(z) -> -q as z->0 when gamma = 0
//
// computed by a consequence of power expansions with improvements near points z=1 and z=\infty
//
// it is assumed that z does not belong to the branch-cut [1,\infty)
//
// Usage:
//  = HeunC(q,alpha,gamma,delta,epsilon,z)
//
// Returned parameters:
// val is the value of the Heun function
// dval is the value of z-derivative of the Heun function
// err is the estimated error
// numb is the number of power series terms needed for the evaluation
//
// Oleg V. Motygin, copyright 2018, license: GNU GPL v3
//
// 26 January 2018
//

#include <algorithm>
#include <complex>
#include <cmath>
#include <vector>
#include <chrono>
#include "Matrix.hpp"
#include "sgn.hpp"

#define M_PI           3.14159265358979323846  /* pi */
#define eps 2.2204e-16

namespace HeunCspace {

	struct HeunCvars {
	std::complex<double>    val;
	std::complex<double>    dval;
	int                numb;
	double             err;
	};
	
	struct HeunCparams {
	std::complex<double> q, alpha, gamma, delta, epsilon;
	};
	
	struct ConnectionVars
	{
	Matrix<std::complex<double>> C10;
	double err;
	int numb;
	};
	
	struct savedataVars
	{
	ConnectionVars Cvars;
	HeunCparams p;
	};
	
	void MakeNan(HeunCvars result){
		// make zero not nan for the moment until I can figure out how c++ nan works
        	result.val = 0;
        	result.dval = 0;
        	result.err = 0;
        	result.numb = 0;
	}

	class HeunC {

   	public:
		// options
		/*------------------------
        	Heun_cont_coef:      for each power expansion of the analytic continuation procedure 
		 	     	the coefficient is the relative (to the radius of convergence) distance from 
 		             	the centre to the calculated point.
	
		Heun_klimit          maximum number of power series' terms.
	
		Heun_optserterms     number of power series' terms considered as in some sense optimal.
	
		Heun_asympt_klimit   maximum number of asymptotic series' terms.
	
		Heun_proxco, Heun_proxcoinf_rel      specifies relative proximity to singular point where special representation is used 
	
		Heun_memlimit	is the maximum number of sets of data (parameters of confluent Heun function and corresponding connection coefficients) which are kept in memory */
		const double Heun_cont_coef = 0.38;
		const int Heun_klimit = 1000;
		const int Heun_optserterms = 40;	
		const int Heun_asympt_klimit = 200;	
		const double Heun_proxco = 0.05, Heun_proxcoinf_rel = 1.0;
		const int Heun_memlimit	= 5;
		
		// data storage vectors (not sure why we need these... )
		std::vector<savedataVars> savedata10, savedata0inf;
		// variables R and N
		bool noRN = true;
		double R = 0;
		double N = 0;
		
		// Constructor 
        	HeunC() {} ;
		
		// main caluclation function
		HeunCvars compute(std::complex<double> alpha, std::complex<double> beta, std::complex<double> gamma, 
                          	std::complex<double> delta, std::complex<double> eta, double z);
	
		// calculation function for the second local solution
        	HeunCvars compute_s(std::complex<double> alpha, std::complex<double> beta, std::complex<double> gamma, 
                          	std::complex<double> delta, std::complex<double> eta, double z);
  	
    	private: 
        	// HeunC0
        	HeunCvars HeunC0(HeunCparams& p, double& z, bool aux=false);	
		HeunCvars HeunC00(HeunCparams& p, double& z, bool aux=false);
		HeunCvars HeunC00gen(HeunCparams& p, double& z);
		HeunCvars HeunC00log(HeunCparams& p, double& z);	
		// HeunCs0
		HeunCvars HeunCs0(HeunCparams& p, double& z);
		HeunCvars HeunCs00(HeunCparams& p, double& z);	
		HeunCvars HeunCs00gamma1(HeunCparams& p, double& z);
		// HeunCfromZ0
		HeunCvars HeunCfromZ0(HeunCparams& p, double& z,double& Z0,std::complex<double>& H0,std::complex<double>& dH0);
		// HeunCconnect
		HeunCvars HeunCconnect(HeunCparams& p, double& z, double& z0,std::complex<double>& H0,std::complex<double>& dH0,double R0=0,bool aux=false);
		// HeunCfaraway
		std::pair<HeunCvars, HeunCvars> HeunCfaraway(HeunCparams& p, double& z);
		// HeunCjoin
		ConnectionVars HeunCjoin0inf(HeunCparams& p,bool aux=false);
		ConnectionVars HeunCjoin10(HeunCparams& p);	
		HeunCvars HeunC1(HeunCparams& p, double& z);	
		HeunCvars HeunCs1(HeunCparams& p, double& z);
		ConnectionVars extrdatfromsav(HeunCparams& p, std::vector<savedataVars>& savedata, bool& consts_known);
		void keepdattosav(savedataVars& s, std::vector<savedataVars>& savedata);	
		// HeunCinf
		HeunCvars HeunCinfA(HeunCparams& p, double& z);
		HeunCvars HeunCinfB(HeunCparams& p, double& z);
		// HeunCutils
		std::complex<double> findcoef4HeunCs(HeunCparams& p);
		void findR();
	};		

#include "HeunC.impl.hpp"
}
	
#endif /* HEUNC_HPP_ */	
