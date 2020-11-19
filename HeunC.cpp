#include "HeunC.hpp"

extern "C" void HeunC_func(double q_real, double q_imag, double alpha_real, double alpha_imag, double gamma_real, double gamma_imag, 
			     double delta_real, double delta_imag, double epsilon_real, double epsilon_imag, double z, 
			     double *val_real, double *val_imag, double *dval_real, double *dval_imag, int *numb, double *err){
	std::complex<double> q(q_real, q_imag);
	std::complex<double> alpha(alpha_real, alpha_imag);
	std::complex<double> gamma(gamma_real, gamma_imag);
	std::complex<double> delta(delta_real, delta_imag);
	std::complex<double> epsilon(epsilon_real, epsilon_imag);
	HeunCspace::HeunC HC;
	auto output = HC.compute(q, alpha, gamma, delta, epsilon, z);
	*val_real = real(output.val);
	*val_imag = imag(output.val);
	*dval_real = real(output.val);
	*dval_imag = imag(output.dval);
	*numb = output.numb;
	*err = output.err;
}
