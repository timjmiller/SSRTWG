#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));
  Type nll=-sum(dnorm(Y,a+b*x,exp(logSigma),true));
  SIMULATE {
    for(int i = 0; i < Y.size(); i++) Y(i) = rnorm(a+b*x(i), exp(logSigma));
    REPORT(Y);
  }
  return nll;
}

