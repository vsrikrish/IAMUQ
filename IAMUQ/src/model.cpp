#include <numeric>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

double update_pop(double P, double Q, NumericVector psi) {
  // pv[0] is the population, pv[1] is the GWP
  double y = Q/P; // compute previous time step per-capita GWP
  
  double pop = P * (1 + (psi[0] * y / (psi[1] + y)) * ((psi[2] - P) / psi[2])); // update population for next time step
  
  return pop;
}

double update_tfp(double A, double alpha, double As) {
  double tfp = A + alpha * A * (1 - (A / As));
  
  return tfp;
}

double update_labor(double P, double pi) {
  double labor = pi * P;
  
  return labor;
}

double update_capital(double Q, double K, double delta, double s) {
  double capital = ((1-delta) * K) + (s * Q);
  
  return capital;
}

double update_gwp(double A, double L, double K, double lambda) {
  double gwp = A * pow(L, lambda) * pow(K, 1-lambda);
  
  return gwp;
}

double update_emis(double Q, NumericVector gamma, NumericVector rho) {
  double phi = std::inner_product(gamma.begin(), gamma.end(), rho.begin(), 0.0);
  double emis = Q * phi;
  
  return emis;
}

// [[Rcpp::export]]
DataFrame model_run(NumericVector yr,
                    double P0,
                    NumericVector psi,
                    double alpha,
                    double A0,
                    double As,
                    double s,
                    double lambda,
                    double delta,
                    double pi,
                    double kappa,
                    NumericMatrix gamma,
                    NumericVector rho,
                    Nullable<NumericVector> init) {
  int n_yr = yr.length();
  
  // initialize storage vectors for each output
  NumericVector P (n_yr);
  NumericVector A (n_yr);
  NumericVector L (n_yr);
  NumericVector K (n_yr);
  NumericVector Q (n_yr);
  NumericVector C (n_yr);
  
  // set initial values
  P[0] = P0;
  L[0] = pi * P[0];

  if (init.isNotNull()) {
    NumericVector x(init);
    Q[0] = x[0];
    C[0] = x[1];
    K[0] = (s / delta) * Q[0];
    A[0] = Q[0] * pow(L[0], -lambda) * pow(K[0], lambda-1);
  } else {
    A[0] = A0;
    K[0] = L[0] * pow(s * A[0] /  delta,  1/lambda);
    Q[0] = update_gwp(A[0], L[0], K[0], lambda);
    C[0] = update_emis(Q[0], gamma(0, _ ), rho);
  }
  
  // loop over years and run model
  for (int i=1; i < n_yr; ++i) {
    P[i] = update_pop(P[i-1], Q[i-1], psi);
    A[i] = update_tfp(A[i-1], alpha, As);
    L[i] = update_labor(P[i], pi);
    K[i] = update_capital(Q[i-1], K[i-1], delta, s);
    Q[i] = update_gwp(A[i], L[i], K[i], lambda);
    C[i] = update_emis(Q[i], gamma(i, _ ), rho);
  }
  
  // form return DataFrame
  DataFrame mout = DataFrame::create(Named("year") = yr,
                                     Named("P") = P,
                                     Named("Q") = Q,
                                     Named("A") = A,
                                     Named("K") = K,
                                     Named("L") = L,
                                     Named("C") = C
                                    );
  
  return mout;
}
