#include "HeunC.hpp"

extern "C" double HeunC_func(double M, double mu, double omega, double a, int l, int m, bool ingoing, bool KS_or_BL, double t, double r){
	KerrBH_Rfunc Rfunc(M, mu, omega, a, l, m);	
	double result = Rfunc.compute(r, t, ingoing, KS_or_BL);
	return result;
}
