#if !defined(HEUNC_HPP_)
#error "This file should only be included through HeunC.hpp"
#endif

#ifndef HEUNC_IMPL_HPP_
#define HEUNC_IMPL_HPP_

// confluent Heun function, the first local solution of the equation
// HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
// at z=0 such that HeunC(0)=1 and
// HeunC'(0)=-q/gamma when gamma is not equal to 0
// HeunC'(z)/log(z) -> -q as z->0 when gamma = 0
//
// computed by a consequence of power expansions with improvements near points z=1 and z=infty
//
// it is assumed that z does not belong to the branch-cut [1,infty)
//
// Usage:
// [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunC(q,alpha,gamma,delta,epsilon,z)
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

inline HeunCvars HeunC::compute(std::complex<double> alpha_, std::complex<double> beta_, std::complex<double> gamma_, 
                          std::complex<double> delta_, std::complex<double> eta_, double z)
{
  HeunCvars result;
 
  HeunCparams p;
  p.q = 0.5*(alpha_*(1 + beta_) - beta_*(1+gamma_) - 2*eta_ - gamma_);
  p.alpha = 0.5*alpha_*(2 + beta_ + gamma_) + delta_;
  p.gamma = beta_ + 1;
  p.delta = gamma_ + 1;
  p.epsilon = alpha_;

  if (z>=1){
    throw std::invalid_argument("HeunC0: z belongs to the branch-cut [1,infty)");
  }
  else {
    findR();
    
    if ( (std::abs(p.epsilon)>1/2)&& (std::abs(p.q) < 2.5) && (std::abs(z)>Heun_proxcoinf_rel*R/(std::abs(eps)+std::abs(p.epsilon))) ) {
      std::pair<HeunCvars, HeunCvars> vars1_vars = HeunCfaraway(p,z);
      result = vars1_vars.second;
    }
    else {
      result = HeunC0(p,z);
    }
    return result; 
  }
}

inline HeunCvars HeunC::compute_s(std::complex<double> alpha_, std::complex<double> beta_, std::complex<double> gamma_, 
                          std::complex<double> delta_, std::complex<double> eta_, double z)
{
  HeunCvars result;
 
  HeunCparams p;
  p.q = 0.5*(alpha_*(1 + beta_) - beta_*(1+gamma_) - 2*eta_ - gamma_);
  p.alpha = 0.5*alpha_*(2 + beta_ + gamma_) + delta_;
  p.gamma = beta_ + 1;
  p.delta = gamma_ + 1;
  p.epsilon = alpha_;

  if (z>=1) {
    throw std::invalid_argument("HeunCfaraway: z belongs to the branch-cut [1,infty)");
  }
  else {
    findR();

    if (( std::abs(p.epsilon)>1/2 )&&( std::abs(p.q)<2.5)&&(std::abs(z)>Heun_proxcoinf_rel*R/(std::abs(eps)+std::abs(p.epsilon)) )) {
      std::pair<HeunCvars, HeunCvars> vars1_vars = HeunCfaraway(p,z);
      result = vars1_vars.second;
    }
    else {
      result = HeunCs0(p,z);
    }
    return result;
  }
}

/* This computes the first local solution around z = 0
*/

// confluent Heun function, the first local solution of the ep.quation
// |z| should not exceed the convergency radius 1
//
// aux is an optional parameter, it only works for p.gamma = 0, -1, -2, ...
// if aux = "yes" then the function is computed as a combination with the second
// solution so as to satisfy
// HeunC00(p.q,p.alpha,p.gamma,p.delta,varp.epsilon;z)=
// exp(-\varp.epsilon z)[HeunC00(p.q-p.epsilon*p.gamma,p.alpha-varp.epsilon(p.gamma+p.delta),p.gamma,p.delta,-varp.epsilon;z)+
// + A * HeunCs00(p.q-p.epsilon*p.gamma,p.alpha-varp.epsilon(p.gamma+p.delta),p.gamma,p.delta,-varp.epsilon;z)]
//
// Returned parameters:
// val is the value of the Heun function
// dval is the value of z-derivative of the Heun function
// err is the estimated error
// numb is the number of power series terms needed for the evaluation
//
// Oleg V. Motygin, copyright 2018, license: GNU GPL v3
//
// 15 February 2018
//

inline HeunCvars HeunC::HeunC0(HeunCparams& p, double& z, bool aux){
  HeunCvars result;
  
  if (z>=1){
    throw std::invalid_argument("HeunC0: z belongs to the branch-cut [1,infty)");
  }
  else {
    bool expgrow = false;
    HeunCparams p1 = p;
    if (aux==false) {
      expgrow = std::real(-p.epsilon*z)>0;
      if (expgrow) {
        p1.q = p.q - p.epsilon*p.gamma;
        p1.alpha = p.alpha - p.epsilon * (p.gamma+p.delta);
        p1.epsilon = -p.epsilon;
      }
    }

    if (std::abs(z)<Heun_cont_coef) {
	result = HeunC00(p1,z,aux);
    }
    else {
      double z0 = Heun_cont_coef*z/std::abs(z);
      HeunCvars result0 = HeunC00(p1,z0,aux);
      HeunCvars result1 = HeunCconnect(p1,z,z0,result0.val,result0.dval,R);
      result.numb = result0.numb + result1.numb;
      result.err = result0.err + result1.err;
      result.val = result1.val;
      result.dval = result1.dval;
    }
    if (expgrow) {
      result.val = result.val * exp(p1.epsilon*z);
      result.dval = (p1.epsilon * result.val + result.dval) * exp(p1.epsilon*z);
      result.err = result.err * std::abs(exp(p1.epsilon*z));
    }
  }
  return result;
}

inline HeunCvars HeunC::HeunC00(HeunCparams& p, double& z, bool aux)
{
	// define the result 
	HeunCvars result;

	// check you are in the convergence range
	if (std::abs(z)>=1){
		throw std::invalid_argument("HeunC00: z is out of the convergence radius = 1");
	}
  	else{
    		bool gamma_is_negative_integer = (std::imag(p.gamma)==0) && std::abs(std::ceil(std::real(p.gamma)-5*eps)+std::abs(p.gamma))<5*eps;
		if (gamma_is_negative_integer) {
	     		result = HeunC00log(p, z);
	      		if (aux) {
	          		std::complex<double> co = findcoef4HeunCs(p);
	          		HeunCvars result_s = HeunCs00(p, z);
		          	result.val = result.val - co * result_s.val;
		          	result.dval = result.dval - co * result_s.dval;
		          	result.err = result.err + std::abs(co) * result_s.err;
		          	result.numb = result.numb + result_s.numb;
			}
		}
    		else{
      		result = HeunC00gen(p, z);
		}
	}
	return result;
}

// confluent Heun function expansion for |z| < 1, gamma is not equal to 0, -1, -2, ...
inline HeunCvars HeunC::HeunC00gen(HeunCparams& p, double& z)
{
	HeunCvars result;

	if (z==0) {
    		result.val = 1;
		result.dval = -p.q/p.gamma;
    		result.err = 0;
		result.numb = 1;
  	}
  	else {
    		// recursion relation variables
    		std::complex<double> ckm2 = 1; 
		std::complex<double> ckm1 = -z*p.q/p.gamma;
    		result.val = ckm2+ckm1; 
		std::complex<double> vm1 = result.val; 
		std::complex<double> vm2;
    		result.dval = -p.q/p.gamma; 
		std::complex<double> dm1 = result.dval; 
		std::complex<double> dm2;
    		std::complex<double> ddval = 0; 
    		int k = 2; 
		std::complex<double> ckm0 = 1;
		// perform recursion
    		while ( (k<=Heun_klimit) && ( (vm2 != vm1) || (dm2 != dm1) || (std::abs(ckm0)>eps) ) ) {
			ckm0 = (ckm1*z*(-p.q+(k-1)*(p.gamma-p.epsilon+p.delta+k-2)) + ckm2*(z*z)*((k-2)*p.epsilon+p.alpha))/(k*(p.gamma+k-1));
		      	result.val = result.val + ckm0; 
			result.dval = dm1 + k*ckm0/z;
      			ddval = ddval + k*(k-1)*ckm0/(z*z);
      			ckm2 = ckm1; 
			ckm1 = ckm0;
      			vm2 = vm1; 
			vm1 = result.val;
      			dm2 = dm1; 
			dm1 = result.dval;
      			k += 1;
		}
    		result.numb = k-1;

		if ( std::isinf(std::abs(result.val)) || std::isinf(std::abs(result.dval)) || std::isnan(std::abs(result.val)) || std::isnan(std::abs(result.dval)) ) {
			throw std::runtime_error("HeunC00: failed convergence of recurrence and summation"); 
      		}
    		else {
			double err1;
			std::complex<double> val2;
			if (p.q-p.alpha*z != 0.0) {
        			val2 = ( z*(z-1)*ddval+(p.gamma*(z-1)+p.delta*z+p.epsilon*z*(z-1))*result.dval ) / (p.q-p.alpha*z);
        			err1 = std::abs(result.val-val2);
      			}
			else {
        			err1 = INFINITY;
      			}
      			if (std::abs(p.q-p.alpha*z)<0.01) {
				double err2;
        			err2 = std::abs(ckm0) * sqrt(result.numb) + eps * result.numb * std::abs(result.val);
        			result.err =  std::min(err1,err2);
      			}
			else {		
        			result.err = err1;
      			}
		}
  	}
	return result;
}

// confluent Heun function, p.gamma = 0, -1, -2, ...
inline HeunCvars HeunC::HeunC00log(HeunCparams& p, double& z) {
	HeunCvars result;
  	
	if (z==0) {
    		result.val = 1; 
    		result.err = 0; 
		result.numb = 1;
    		if (std::abs(p.gamma)<eps) {
     			result.dval = INFINITY;
    		}
		else {
      		result.dval = -p.q/p.gamma;
    		}
  	}
	else {
    		int N = std::round(1-std::real(p.gamma));
  		
		// recursion relation variables
    		std::complex<double> L1 = 1, dL1 = 0, ddL1 = 0;
		std::complex<double> L2 = 0, dL2 = 0, ddL2 = 0;
		std::complex<double> L3, dL3, ddL3; 
		std::complex<double> ddval;
		std::complex<double> ckm0 = 1, ckm1 = 1, ckm2 = 0;
		std::complex<double> dm1, dm2, skm0, skm1, skm2 = 0;
		std::complex<double> dsm1, dsm2;

		for(int k=1; k<N; k++) {
			ckm0 = (ckm1*z*(-p.q+(k-1)*(p.gamma-p.epsilon+p.delta+k-2)) + ckm2*(z*z)*((k-2)*p.epsilon+p.alpha))/(k*(p.gamma+k-1));
	      		L1 = L1+ckm0; 
			dL1 = dL1+k*ckm0/z; 
			ddL1 = ddL1+k*(k-1)*ckm0/(z*z);
	      		ckm2 = ckm1; 
			ckm1 = ckm0; 
	    	}
	       
                std::complex<double> sN = (ckm1*z*(p.q+p.gamma*(p.delta-p.epsilon-1)) + ckm2*(z*z)*(p.epsilon*(p.gamma+1)-p.alpha))/(p.gamma-1);

		dm1 = dL2; 
		//dm2 = nan;
	    	ckm1 = 0; 
		ckm2 = ckm0; 
	    	L3 = sN; 
		skm1 = sN;
	    	dL3 = N*sN/z; 
		ddL3 = N*(N-1)*sN/(z*z); 
	    	dsm1 = dL3; 
		//dsm2 = nan; 
		//skm0 = nan;
	    
	    	int k = N+1;

	    	while ( (k<=Heun_klimit) && ( (dsm2!=dsm1) || (std::abs(skm0)>eps) || (dm2!=dm1) || (std::abs(ckm0)>eps) ) ) {
	      		skm0 = (skm1*z*(-p.q+(k-1)*(p.gamma-p.epsilon+p.delta+k-2)) + skm2*(z*z)*((k-2)*p.epsilon+p.alpha))/(k*(p.gamma+k-1));
	      		ckm0 = (ckm1*z*(-p.q+(k-1)*(p.gamma-p.epsilon+p.delta+k-2)) + ckm2*(z*z)*((k-2)*p.epsilon+p.alpha))/(k*(p.gamma+k-1)) +
				(-skm0*(p.gamma+2*k-1)+skm1*z*(p.gamma-p.epsilon+p.delta+2*k-3)+skm2*(z*z)*p.epsilon)/(k*(p.gamma+k-1));
	      		L2 = L2+ckm0; 
			dL2 = dm1+k*ckm0/z; 
			ddL2 = ddL2+k*(k-1)*ckm0/(z*z);
			ckm2 = ckm1; 
			ckm1 = ckm0;
			dm2 = dm1; 
			dm1 = dL2;

	      		L3 = L3+skm0; dL3 = dsm1+k*skm0/z; ddL3 = ddL3+k*(k-1)*skm0/(z*z);
	      		skm2 = skm1; skm1 = skm0;
	      		dsm2 = dsm1; dsm1 = dL3;
	      		k++;
		}
	    	result.numb = k-1;
	   	result.val = L1 + L2 + std::log(z) * L3;
	    	result.dval = dL1 + dL2 + std::log(z) * dL3 + L3/z;
	    	ddval = ddL1 + ddL2 - L3/(z*z) + 2*dL3/z + std::log(z) * ddL3;

	    	if ( ( std::isinf(std::abs(result.val)) || std::isinf(std::abs(result.dval)) ) || ( std::isnan(std::abs(result.val)) || std::isnan(std::abs(result.dval)) ) ) {
                        throw std::runtime_error("HeunC00log: failed convergence of recurrence and summation"); 
                }
	    	else {
			std::complex<double> val2, val3;
	      		double err1, err2;
			if (p.q-p.alpha*z != 0.0) {	    
	       			val2 = ( z*(z-1)*ddval+(p.gamma*(z-1)+p.delta*z+p.epsilon*z*(z-1))*result.dval ) / (p.q-p.alpha*z);
	        		val3 = ((dL3*p.epsilon+ddL3)*(z*z)*std::log(z)+(dL3*(p.gamma-p.epsilon+p.delta)-ddL3)*z*std::log(z)-dL3*p.gamma*std::log(z)+
	          			(p.epsilon*(dL2+dL1)+ddL2+ddL1)*(z*z)+((dL1+dL2)*(p.gamma-p.epsilon+p.delta)+L3*p.epsilon-ddL2-ddL1+2*dL3)*z+
	          			L3*(1-p.gamma)/z-(dL1+dL2)*p.gamma+L3*(p.gamma+p.delta-p.epsilon)-2*dL3-L3) / (p.q-p.alpha*z);
	      
	        		err1 = std::min(std::abs(result.val-val2),std::abs(result.val-val3));
			}
	      		else {
	        		err1 = INFINITY;
	      		}
	      		if ((std::abs(p.q-p.alpha*z)<0.01)||(err1<eps)) {
	        		err2 = std::abs(L1)*eps*N + std::abs(ckm0)*sqrt(result.numb-N+1) + std::abs(L2)*eps*(result.numb-N+1) +
	               			std::abs(std::log(z)) * ( std::abs(skm0)*sqrt(result.numb-N+1) + std::abs(L3)*eps*(result.numb-N+1) );
	        		result.err =  std::min(err1,err2);
	      		}
			else {
	        		result.err = err1;
	      		}
		}
	}
	return result;
}

// confluent Heun function,
// the second local solution
// with branch-cut (1,+\infinity)
// HeunCs(z) = z^(1-p.gamma)*h(z), where h(0)=1, 
// h'(0)=(-q+(1-p.gamma)*(p.delta-p.epsilon))/(2-p.gamma) for p.gamma not equal to 1, 2
// h'(z)/log(z) -> -q+(1-p.gamma)*(p.delta-p.epsilon) as z\to0 for p.gamma=2
// and
// HeunCs(z) \sim log(z) - q * z * log(z) +   as z\to0 for p.gamma=1
//
//
// Usage:
// [val,result.dval,err,numb,wrnmsg] = HeunCs00(q,p.alpha,p.gamma,p.delta,p.epsilon,z)
//
// Returned parameters:
// val is the value of the confluent Heun function
// dval is the value of z-derivative of the confluent Heun function
// err is the estimatedresult.error
// numb is the total result.number of power series terms needed for the evaluation
//
// Oleg V. Motygin, copyright 2018, license: GNU GPL v3
//
// 09 January 2018
//

// the second local solution at z=0 (see HeunCs00)
//
// computed by a consequence of power expansions
inline HeunCvars HeunC::HeunCs0(HeunCparams& p, double& z){
  
  HeunCvars result;

  if (z>=1){
    throw std::invalid_argument("HeunC0: z belongs to the branch-cut [1,infty)");
  }
  else {
    HeunCparams p1 = p;
    bool expgrow = std::real(-p.epsilon*z)>0;
    if (expgrow) {
      p1.q = p.q - p.epsilon * p.gamma;
      p1.alpha = p.alpha - p.epsilon * (p.gamma+p.delta);
      p1.epsilon = -p.epsilon;
    } 
    
    if (std::abs(z)<Heun_cont_coef){
      result = HeunCs00(p1,z);
    }
    else {
      double z0 = Heun_cont_coef*z/std::abs(z);
      HeunCvars result0 = HeunCs00(p1,z0);
      HeunCvars result1 = HeunCconnect(p1,z,z0,result0.val,result0.dval,R);
      result.numb = result0.numb + result1.numb;
      result.err = result0.err + result1.err;
    }
    if (expgrow) {
      result.val = result.val * exp(p.epsilon*z);
      result.dval = (p.epsilon*result.val + result.dval) * exp(p.epsilon*z);
      result.err = result.err * std::abs(exp(p.epsilon*z));
    }
    return result;
  }
}

// solution at z ~ 0
// |z| should not exceed the convergency radius 1
inline HeunCvars HeunC::HeunCs00(HeunCparams& p, double& z)
{
  HeunCvars result;
  if (std::abs(z)>=1){
        throw std::invalid_argument("HeunCs00: z is outside the |z|<1 radius of convergence");
  }

  else {
    if (std::abs(p.gamma-1)<eps){
      if ( z==0 ){
       result.val= INFINITY; result.dval = INFINITY;
       result.err = INFINITY; result.numb = 1;
      }
      else {
	result = HeunCs00gamma1(p,z);
      }
    } 
    else {
      HeunCparams p1;
      p1 = p;
      p1.q = p.q + (p.gamma-1)*(p.delta-p.epsilon);
      p1.alpha = p.alpha+p.epsilon*(1-p.gamma);
      p1.gamma = 2-p.gamma;
      std::cout << " using HeunC00 z = " << z << std::endl;
      HeunCvars H0 = HeunC00(p1,z);
      result.val= std::pow(z,(1-p.gamma)*H0.val);
      result.dval = (1-p.gamma)*std::pow(z,(-p.gamma)*H0.val) + std::pow(z,(1-p.gamma)*H0.dval);
      if ( std::isinf(std::abs(result.val)) || std::isinf(std::abs(result.dval)) ){
         result.err = INFINITY;
      }
      else {
         result.err = std::abs(std::pow(z,(1-p.gamma))*H0.err);
      }
    }
    return result; 
  }
}

// confluent Heun function, second local solution at z=0, p.gamma = 1
//
inline HeunCvars HeunC::HeunCs00gamma1(HeunCparams& p,double& z)
{  
  HeunCvars result;

  // declare iteration variables
  std::complex<double> L1 = 0, dL1 = 0, ddL1 = 0, dm1 = 0, dm2, ckm0, ckm1 = 0, ckm2 = 0; 
  std::complex<double> L2 = 1, dL2 = 0, ddL2 = 0, skm2 = 0, skm1 = 1, dsm1 = 0, dsm2, skm0;
  std::complex<double> ddval;
  int k = 1;

  while ( (k<=Heun_klimit) && ( (dsm2!=dsm1) || (std::abs(skm0)>eps) || (dm2!=dm1) || (std::abs(ckm0)>eps) ) ){

    skm0 = (skm1*z*(-p.q+(k-1)*(-p.epsilon+p.delta+k-1)) + skm2*(z*z)*((k-2)*p.epsilon+p.alpha))/(k*k);

    ckm0 = (ckm1*z*(-p.q+(k-1)*(-p.epsilon+p.delta+k-1)) + ckm2*(z*z)*((k-2)*p.epsilon+p.alpha))/(k*k)+
           (skm0*2 + skm1*z*(p.epsilon/k-p.delta/k-2+2/k))/k + skm2*(z*z)*p.epsilon/(k*k);

    L1 = L1+ckm0; dL1 = dm1+k*ckm0/z; ddL1 = ddL1+k*(k-1)*ckm0/(z*z);
    ckm2 = ckm1; ckm1 = ckm0;
    dm2 = dm1; dm1 = dL1;
    
    L2 = L2+skm0; dL2 = dsm1+k*skm0/z; ddL2 = ddL2+k*(k-1)*skm0/(z*z);
    skm2 = skm1; skm1 = skm0;
    dsm2 = dsm1; dsm1 = dL2;
    k += 1;
  }
  result.numb = k-1;
  result.val = L1 + std::log(z) * L2;
  result.dval = dL1 + std::log(z) * dL2 + L2/z;
  ddval = ddL1 - L2/(z*z) + 2*dL2/z + std::log(z) * ddL2;

  if ( ( std::isinf(std::abs(result.val)) || std::isinf(std::abs(result.dval)) ) || ( std::isnan(std::abs(result.val)) || std::isnan(std::abs(result.dval)) ) ){
     throw std::runtime_error("HeunCs00gamma1: failed convergence of recurrence and summation"); 
  }
  else {
    std::complex<double> val2, val3;
    double err1;
    if (p.q-p.alpha*z!= 0.0){
  
      val2 = ( z*(z-1)*ddval+(z-1+p.delta*z+p.epsilon*z*(z-1))*result.dval) / (p.q-p.alpha*z);
    
      val3 = ((dL2*p.epsilon+ddL2)*(z*z)*log(z)+(dL2*(-p.epsilon+p.delta+1)-ddL2)*z*log(z)-dL2*log(z)+ 
             (dL1*p.epsilon+ddL1)*(z*z)+(dL1*(-p.epsilon+p.delta+1)+L2*p.epsilon-ddL1+2*dL2)*z-L2*p.epsilon+ 
             L2*p.delta-2*dL2-dL1) / (p.q-p.alpha*z);
             
      err1 = std::min(std::abs(result.val-val2),std::abs(result.val-val3));   
    }
    else {
     err1 = INFINITY;
    }
    if ((std::abs(p.q-p.alpha*z)<0.01)||(err1<eps)){
     double err2 = std::abs(ckm0)*sqrt(result.numb) + std::abs(L1)*eps*result.numb + 
             std::abs(log(z)) * ( std::abs(skm0)*sqrt(result.numb) + std::abs(L2)*eps*result.numb );           
     result.err = err2+std::min(err1,err2);
    }
    else {
     result.err = err1;
    }
    return result;
  }
}


// confluent Heun function, a solution of the ep.quation
// HeunC""(z)+(p.gamma/z+p.delta/(z-1)+p.epsilon)*HeunC"(z)+(p.alpha*z-p.q)/(z*(z-1))*HeunC(z) = 0
// computed at z by power series about Z0 for the given values H(Z0)=H0, H"(Z0)=dH0 
//
// it is assumed that z, Z0 are not equal to 0, 1 and |z-Z0| < std::min{|Z0|,|Z0-1|}
//
// Usage:
// [val,dval,err,numb] = HeunCfromZ0(p.q,p.alpha,p.gamma,p.delta,p.epsilon,z,Z0,H0,dH0)
//
// Returned parameters:
// val is the value of the confluent Heun function at point z
// dval is the value of z-derivative of the Heun function at point z
// err is the estimated error
// numb is the number of power series terms needed for the evaluation
//
// Oleg V. Motygin, copyright 2018, license: GNU GPL v3
//
// 09 January 2018
//
inline HeunCvars HeunC::HeunCfromZ0(HeunCparams& p,double& z,double& Z0,std::complex<double>& H0,std::complex<double>& dH0)
{
  HeunCvars result;

  R = std::min(std::abs(Z0),std::abs(Z0-1.0));
  
  if (std::abs(z-Z0)>=R) {
    throw std::invalid_argument("HeunCfromZ0: z is out of the convergence radius"); 
  }
  else if ((std::abs(z-1)<eps) || (std::abs(Z0-1)<eps)) {
   throw std::invalid_argument("HeunCfromZ0: z or Z0 is too close to the singular points"); 
  }
  else if (z==Z0) {
   result.val= H0; result.dval = dH0; 
   result.err= 0; result.numb = 0;
  } 
  else {
    double zeta = z-Z0;
    // iteration variables
    std::complex<double> ckm0, ckm1, ckm2, ckm3;
    std::complex<double> dm1, dm2, vm1, vm2;  
    std::complex<double> ddval;
    ckm3 = H0; 
    ckm2 = dH0*zeta; 

    // initialise with long recursion relation
    ckm1 = (ckm2*zeta*(2-1)*(p.epsilon*(Z0*Z0)+(p.gamma-p.epsilon+p.delta)*Z0-p.gamma)+
        ckm3*std::pow(zeta,2)*((p.alpha)*Z0-p.q)) / (Z0*(Z0-1)*(1-2)*2);

    result.val= ckm3 + ckm2 + ckm1; 
    vm1 = result.val;
    dm2 = dH0; dm1 = dH0 + 2*ckm1/zeta; 
    result.dval = dm1;
    ddval = 2*ckm1/std::pow(zeta,2); 
  
    int k = 3; 
    ckm0 = 1;
    
    while ( (k<=Heun_klimit) && ( ( vm2!=vm1 ) || ( dm2!=dm1 ) || (std::abs(ckm0)>eps) ) ) {
      // long recursion relation
      ckm0 = (ckm1*zeta*(k-1)*(p.epsilon*(Z0*Z0)+(p.gamma-p.epsilon+p.delta+2*(k-2))*Z0-p.gamma-k+2)+ 
             ckm2*std::pow(zeta,2)*((2*(k-2)*p.epsilon+p.alpha)*Z0-p.q+(k-2)*(p.gamma-p.epsilon+p.delta+k-3))+ 
             ckm3*std::pow(zeta,3)*((k-3)*p.epsilon+p.alpha)) / (Z0*(Z0-1)*(1-k)*k);
      result.val += ckm0; 
      result.dval = dm1 + k*ckm0/zeta;
      ddval += k*(k-1)*ckm0/std::pow(zeta,2);
      ckm3 = ckm2; 
      ckm2 = ckm1; 
      ckm1 = ckm0;
      vm2 = vm1; vm1 = result.val;
      dm2 = dm1; dm1 = result.dval;
      k++;
    }
    
    result.numb = k-1;

    if ( std::isinf(std::abs(result.val)) || std::isinf(std::abs(result.dval)) || std::isnan(std::abs(result.val)) || std::isnan(std::abs(result.dval)) ) {
      throw std::runtime_error("HeunCfromZ0: failed convergence of recurrence and summation"); 
    }
    else {
      std::complex<double> val2;
      double err1, err2;
      if (p.q-p.alpha*z != 0.0) {
        val2 = ( z*(z-1)*ddval+(p.gamma*(z-1)+p.delta*z+p.epsilon*z*(z-1))*result.dval ) / (p.q-p.alpha*z);
        err1 = std::abs(result.val-val2);
      }
      else {
        err1 = INFINITY;
      }
      if (std::abs(p.q-p.alpha*z)<0.01) {
        err2 = std::abs(ckm0) * sqrt(result.numb) + std::abs(result.val) * eps * result.numb;
        result.err =  std::min(err1, err2);
      }
      else {
       result.err = err1;
      }
    }
  }
  return result; 
}

// HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
// by analytic continuation from point z0, where HeunC(z0) = H0, HeunC'(z0) = dH0,
// to another point z, along the line [z0,z], using a consequence of power expansions
//
// Assumptions:
// z0, z are not 0 or 1
//
// R0 is an optional parameter, step size's guess
//
// R is the size of the last used step
//
// Oleg V. Motygin, copyright 2018, license: GNU GPL v3
//
// 09 January 2018
//

inline HeunCvars HeunC::HeunCconnect(HeunCparams& p, double& z, double& z0,std::complex<double>& H0,std::complex<double>& dH0, double R0, bool aux)
{
  HeunCvars result;

  if ((z==0)||(z==1)||(z0==0)||(z0==1)){
     throw std::invalid_argument("HeunCconnect: assumed that z, z0 are not equal to 0, 1");
  }
  else {
    bool insearch = true;
    double Rmax;    

    if (R0>0) {
	Rmax = R0;
    } 
    else {
      Rmax = 12*Heun_cont_coef/(1+std::abs(p.epsilon));
    }

    R = std::min(Rmax,std::min(std::abs(z0),std::abs(z0-1))*Heun_cont_coef);

    // first set of iteration variables
    bool Rtuned = false;
    int iter = 1;
    double z1;       // step end point
    int positivity = 2*(z >= z0) - 1;

    while (Rtuned==false) {
      if (std::abs(z-z0) <= R) {
        z1 = z;
      }
      else {
        z1 = z0 + R * positivity;
      }
      
      result = HeunCfromZ0(p,z1,z0,H0,dH0);

      Rtuned = (result.err < 5*eps) && (result.numb < Heun_optserterms) || (iter>5) || (result.numb<=8);
      
      if (Rtuned==false) {
        R = R / std::max(result.err/(5*eps), static_cast<double>(result.numb)/Heun_optserterms);
      }

      insearch = !(Rtuned && (z==z1));
      iter = iter+1; 
    }

    //second set of iteration variables
    double errsum;
    int numbsum;

    z0 = z1;
    errsum = result.err; 
    numbsum = result.numb; 
    HeunCvars result_; // iteration version of result
    H0 = result.val; dH0 = result.dval;
    
    while (insearch) {
      R = std::min(R, std::min(std::abs(z0),std::abs(z0-1))*Heun_cont_coef);
      if (std::abs(z-z0) <= R) {
        z1 = z; insearch = false;
      }
      else {
        z1 = z0 + R * positivity;
      }
      result_ = HeunCfromZ0(p,z1,z0,H0,dH0);
      H0 = result_.val;
      dH0 = result_.dval;

      errsum += result_.err;
      numbsum += result_.numb;

      if (insearch) {
        R = Heun_optserterms * R / (result_.numb + eps);
      }
      
      z0 = z1;        
    }
    result.numb = numbsum; 
    result.val = H0; result.dval = dH0; result.err = errsum;
  }
  return result;
}

/*
Depends on:

HeunCconnect
HeunCjoin0inf
HeunCinfA
HeunCinfB

*/
// confluent Heun function, a solution of the equation
// HeunC""(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC"(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
//
// computation for  sufficiently large |z|, by analytic continuation from infinity
//
// computes both the first at z=0 local solution (see HeunC00) and the second at z=0 local solution (see HeunCs0)
//
// It is assumed that epsilon \neq 0 !
//
// Usage:
// [val1,dval1,result1.err,result2.val,result2.dval,result2.err,numb,,val1woexp,dval1woexp,result1.errwoexp,result2.valwoexp,result2.dvalwoexp,result2.errwoexp] = HeunCfaraway(p,z)
//
// Returned parameters:
// val1 is the value of the Heun function, growing from the first local solution at z=0
// dval1 is the value of z-derivative of the Heun function
// result1.err is the estimated error
// result2.val is the value of the Heun function, growing from the second local solution at z=0
// result2.dval is the value of z-derivative of the Heun function
// result2.err is the estimated error
// numb is the number of power series terms needed for the evaluation
//
// val1woexp = val1 * exp(epsilon*z)
// dval1woexp = dval1 * exp(epsilon*z)
// result1.errwoexp = result1.err * std::abs(exp(epsilon*z))
//
// result2.valwoexp = result2.val * exp(epsilon*z)
// result2.dvalwoexp = result2.dval * exp(epsilon*z)
// result2.errwoexp = result2.err * std::abs(exp(epsilon*z))
//
// Oleg V. Motygin, copyright 2018, license: GNU GPL v3
//
// 26 March 2018
//

//

// for large |z|
inline std::pair<HeunCvars, HeunCvars> HeunC::HeunCfaraway(HeunCparams& p, double& z)
{
  HeunCvars result1, result2;

  if (z>=1) {
    throw std::invalid_argument("HeunCfaraway: z belongs to the branch-cut [1,infty)");
  }
  else {
    HeunCparams pB = p;
    pB.q = p.q-p.epsilon*p.gamma;
    pB.alpha = p.alpha-p.epsilon*(p.gamma+p.delta);
    pB.epsilon = -p.epsilon;

    ConnectionVars CA = HeunCjoin0inf(p);
    ConnectionVars CB = HeunCjoin0inf(pB);

    findR();

    double infpt = -std::max(1.0,R/(std::abs(eps)+std::abs(p.epsilon)));

    HeunCvars varsA, varsB;

    if (std::abs(z)>std::abs(infpt)) {
      
      HeunCvars varsA = HeunCinfA(p,z);
      HeunCvars varsB = HeunCinfA(pB,z);
    
      result1.numb = CA.numb + CB.numb + varsA.numb + varsB.numb;
      result1.err = CA.err + CB.err + varsA.err + varsB.err;
      result2.numb = result1.numb;    
      result2.err = result1.err;
    } 
    else { 
      HeunCvars varsinfA = HeunCinfA(p,z);
      varsA = HeunCconnect(p,z,infpt,varsinfA.val,varsinfA.dval, R);

      HeunCvars varsinfB = HeunCinfA(pB,z);
      varsB = HeunCconnect(pB,z,infpt,varsinfA.val,varsinfA.dval, R);
    
      result1.numb = varsinfA.numb + varsinfB.numb + varsA.numb + varsB.numb;
      result1.err = varsinfA.err + varsinfB.err + varsA.err + varsB.err;
      result2.err = result1.err;
    }
    Matrix<std::complex<double>> m;
    m.init(CA.C10(0, 0),CA.C10(1, 0),CB.C10(0, 0),CB.C10(1, 0));
    double co = m.cond();

    result1.val = m(0, 0) * varsA.val + m(0, 1) * exp(-p.epsilon*z) * varsB.val;
    result1.dval = m(0, 0) * varsA.dval + m(0, 1) * exp(-p.epsilon*z) * (-p.epsilon * varsB.val + varsB.dval);
    result1.err = co * (std::abs(m(0, 0)) * varsA.err + std::abs(m(0, 1)) * std::abs(exp(-p.epsilon*z)) * varsB.err);
    
    result2.val = m(1, 0) * varsA.val + m(1, 1) * exp(-p.epsilon*z) * varsB.val;
    result2.dval = m(1, 0) * varsA.dval + m(1, 1) * exp(-p.epsilon*z) * (-p.epsilon * varsB.val + varsB.dval);
    result2.err = co * (std::abs(m(1, 0)) * varsA.err + std::abs(m(1, 1)) * std::abs(exp(-p.epsilon*z)) * varsB.err);
    
    std::pair<HeunCvars, HeunCvars> output(result1, result2);
    return output;
  }
}

// for confluent Heun function, a solution of the equation
// HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
//
// HeunCjoin0inf finds connection coefficients C0, Cs, such that
// C0 * HeunC00(z) + Cs * HeunCs00(z) analytically continues to
// the first, power solution at infinity \exp(i\theta) infty
// (see HeunCinfA)
//
// Returned parameters:
// C10 is the matrix of connection coefficients
// err is the estimated error
// numb is the number of power series terms needed for the evaluation
//
// Oleg V. Motygin, copyright 2018, license: GNU GPL v3
//
// 15 March 2018
//

inline ConnectionVars HeunC::HeunCjoin0inf(HeunCparams& p,bool aux)
{
  ConnectionVars result;
  bool consts_known = false;

  result = extrdatfromsav(p, savedata0inf, consts_known);

  if (consts_known) {
    result.numb = 0;
  }
  else {
    HeunCvars varsinf, varsJinf, varsJ0, varsJs;

    findR();

    double R0, infpt, joinpt;

    R0 = R/(std::abs(eps)+std::abs(p.epsilon));

    infpt = -2 * R0;
    joinpt = -std::min(1.0,R0);

    varsinf = HeunCinfA(p,infpt);       // value at "infinity" point
    varsJinf = HeunCconnect(p,joinpt,infpt,varsinf.val,varsinf.dval); //

    varsJ0 = HeunC0(p,joinpt,aux); // first solution at join point near zero 
    varsJs = HeunCs0(p,joinpt);         // second solution at join point near zero

    result.err = varsinf.err + varsJinf.err + varsJ0.err + varsJs.err;
    result.numb = varsinf.numb + varsJinf.numb + varsJ0.numb + varsJs.numb;

    Matrix<std::complex<double>> m, b;
    m.init(varsJ0.val,varsJs.val, varsJ0.dval,varsJs.dval);
    b.init(varsJinf.val, 0, varsJinf.dval, 0);
    result.C10 = (m.inverse().dot(b));

    savedataVars s;
    s.p = p;
    s.Cvars = result;
    keepdattosav(s, savedata0inf);
  }
  return result;
}

inline ConnectionVars HeunC::extrdatfromsav(HeunCparams& p, std::vector<savedataVars>& savedata, bool& consts_known){
  ConnectionVars result;
  result.numb = 0; 
  savedataVars s;

  if (savedata.size() !=0) {
    for(int k=1; k < savedata.size(); k++){
        s = savedata[k];
        if (s.p.q == p.q && s.p.alpha == p.alpha && s.p.gamma == p.gamma && s.p.delta == p.delta && s.p.epsilon == p.epsilon) {
        	result = s.Cvars;        
        	consts_known = "true";
        	break;
      	}
    }
  }
  return result;
}

inline void HeunC::keepdattosav(savedataVars& s, std::vector<savedataVars>& savedata)
{
  if (savedata.size() > Heun_memlimit)
  {
    savedata.erase(savedata.begin()); // remove first element
  }
  savedata.push_back(s); // store new element  
}

// confluent Heun function
// asymptotic expansion at z=infinity
// the first, power solution
//
// Usage:
// [val,dval,err,numb] = HeunCinfA(p.q,p.alpha,p.gamma,p.delta,p.epsilon,z)
//
// Returned parameters:
// val is the value of the Heun function
// dval is the value of z-derivative of the Heun function
// err is the estimated error
// numb is the number of the summed series terms
//
// Oleg V. Motygin, copyright 2017-2018, license: GNU GPL v3
//
// 20 December 2017
//

inline HeunCvars HeunC::HeunCinfA(HeunCparams& p, double& z)
{
  HeunCvars result;

  result.val = 1; 
  result.dval = 0;
  result.err = 0; 
  result.numb = 1;
  
  // set up iteration variables
  std::complex<double> cnm0, cnm1, cnm2 = 1, cnm3; 
  std::complex<double> dnm0, dnm1, dnm2 = 0, dnm3;
  std::complex<double> vm0, vm1, vm2, vm3; 
  std::complex<double> dvm0, dvm1, dvm2, dvm3;

  cnm1 = cnm2/(z*p.epsilon)*(1+(-p.q+p.alpha/p.epsilon*(2-p.gamma-p.delta-1+p.alpha/p.epsilon)+p.alpha-1));
  dnm1 = -cnm1/z;

  result.val = cnm2 + cnm1; 
  result.dval = dnm1;
  
  vm0 = result.val;
  dvm0 = result.dval;  

  result.numb = 2; 
  double small = sqrt(eps); 
  
  bool growcn = false, growdn = false, valstab = false, dvalstab = false;
  int n;  

  while ( (result.numb<=Heun_asympt_klimit) && ((std::abs(cnm3)>small) || !(growcn||valstab) || !(growdn||dvalstab)) ){
    n = result.numb;
    cnm0 = cnm1*n/(z*p.epsilon)*(1+(-p.q+p.alpha/p.epsilon*(2*n-p.gamma-p.delta-1+p.alpha/p.epsilon)+
           (p.gamma-p.epsilon+p.delta+1)*(1-n)+p.alpha-1)/(n*n)) + cnm2/((z*z)*p.epsilon)*
           ((n-2+p.alpha/p.epsilon)*(p.gamma-n+1-p.alpha/p.epsilon))/n;

    dnm0 = -result.numb*cnm0/z;
    result.val += cnm0; 
    result.dval = dnm0;
    result.err = std::abs(cnm2);
    result.numb += 1;
    
    growcn = growcn || ((std::abs(cnm0)>std::abs(cnm1))&&(std::abs(cnm1)>std::abs(cnm2))&&(std::abs(cnm2)>std::abs(cnm3)));
    valstab = valstab || ((vm3==vm2)&&(vm2==vm1)&&(vm1==result.val));

    growdn = growdn || ((std::abs(dnm0)>std::abs(dnm1))&&(std::abs(dnm1)>std::abs(dnm2))&&(std::abs(dnm2)>std::abs(dnm3)));
    dvalstab = dvalstab || ((dvm3==dvm2)&&(dvm2==dvm1)&&(dvm1==result.dval));
    
    if ((std::abs(cnm2)>small) || !(growcn||valstab)){
      cnm3 = cnm2; cnm2 = cnm1; cnm1 = cnm0;
      vm3 = vm2; vm2 = vm1; vm1 = vm0; vm0 = result.val;
    }
  
    if ((std::abs(cnm2)>small) || !(growdn||dvalstab)){
      dnm3 = dnm2; dnm2 = dnm1; dnm1 = dnm0;
      dvm3 = dvm2; dvm2 = dvm1; dvm1 = dvm0; dvm0 = result.dval;
    }
    
  }  
  
  result.val = std::pow((-z),(-p.alpha/p.epsilon)) * vm3;
  result.dval = std::pow((-z),(-p.alpha/p.epsilon)) * (dvm3-p.alpha/p.epsilon*vm3/z);
  result.err = std::abs(std::pow(z,(-p.alpha/p.epsilon))) * result.err;
  
  if ( std::isinf(std::abs(result.val)) || std::isinf(std::abs(result.dval)) || std::isnan(std::abs(result.val)) || std::isnan(std::abs(result.dval)) ){
    throw std::runtime_error("HeunCinfA: failed convergence of recurrence and summation"); 
  }
  return result;
}

// asymptotic expansion at z=infinity,
// the second solution, including exponential factor

inline HeunCvars HeunC::HeunCinfB(HeunCparams& p, double& z)
{  
  HeunCvars result0, result;
  HeunCparams p0 = p;
  p0.q = p.q - p.epsilon*p.gamma;
  p0.alpha = p.alpha - p.epsilon*(p.gamma + p.delta);
  p0.epsilon = -p.epsilon;

  result0 = HeunCinfA(p0,z);
  result = result0;  

  result.val = exp(-p.epsilon*z) * result0.val;
  result.dval = exp(-p.epsilon*z) * (-p.epsilon * result0.val + result0.dval);
  result.err = std::abs(exp(-p.epsilon*z)) * result0.err;
  
  return result;
}    

/* 

Some extra utiliy functions used in the HeunC code

*/

// function find coefficient for the second HeunC solution
inline std::complex<double> HeunC::findcoef4HeunCs(HeunCparams& p){
        int n = std::round(1-std::real(p.gamma));  
        std::complex<double> ckm1 = 1; 
        std::complex<double> ckm2 = 0;
        std::complex<double> ckm0;
        std::complex<double> co = std::pow(p.epsilon,n)/std::tgamma(n+1);
 
        for(int k=1; k < n; k++) {
                ckm0 = (ckm1*(-p.q+(k-1)*(p.gamma-p.epsilon+p.delta+k-2)) + ckm2*((k-2)*p.epsilon+p.alpha))/(k*(p.gamma+k-1));
                co = co + ckm0 * std::pow(p.epsilon,(n-k))/std::tgamma(n-k+1); // tgamma(n+1) = factorial(n)
                ckm2 = ckm1; 
                ckm1 = ckm0; 
        }
        return co;
}

// in HeunC*, HeunCs*, it is assumed that for |z|>2*R, asymptotic series 
// of type \sum_{n=0}^{\infty} b_n n!/z^n can be computed as superasymptotic
// with accuracy better than machine epsilon

inline void HeunC::findR()
{
  if (noRN) {
    double R0;
    double logeps = -36.043674;
    R = -logeps;
    double logfact = 0;
    int n = 1;
  
    while (true)
    {
      n += 1;
      logfact += std::log(static_cast<double>(n));
      R0 = R;
      R = (logfact-logeps)/n;
      if (R > R0){
        break;
      }
    }
    noRN = false;
    N = n-1;
    //R = static_cast<double>(N)*std::exp(-logeps/N)/std::exp(1);
    R = std::exp((logfact - std::log(static_cast<double>(n)) -logeps)/N);
    //R = std::pow((fact/n/eps),(1.0/N));
  }
}

#endif /* HEUNC_IMPL_HPP_ */
