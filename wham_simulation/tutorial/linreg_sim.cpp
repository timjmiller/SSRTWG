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
  vector<Type> Yhat = a + b*x;
  REPORT(Yhat);
  Type nll=-sum(dnorm(Y,Yhat,exp(logSigma),true));
  SIMULATE {
    for(int i = 0; i < Y.size(); i++) Y(i) = rnorm(Yhat(i), exp(logSigma));
    REPORT(Y);
  }
  return nll;
}

