#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix sym_bind_mat(List M, IntegerMatrix H) {
  NumericMatrix M0 = M[0];
  int nr_m = M0.nrow();
  int nc_m = M0.ncol();
  int nr_h = H.nrow();
  int nc_h = H.ncol();
  
  NumericMatrix S = no_init_matrix( nr_m * nr_h, nc_m * nc_h);
  for (int i=0; i < nr_h; ++i) {
    for (int j=0; j < nc_h; ++j) {
      int lind = H(i, j);
      NumericMatrix mat = M[lind];
      NumericMatrix::Sub s = S(Range(nr_m * i, nr_m * (i+1) - 1), Range(nc_m * j, nc_m * (j+1) - 1));
      for (int k=0; k < mat.nrow(); ++k) {
        for (int l=0; l < mat.ncol(); ++l) {
          if (i >= j) {
            s(k, l) = mat(k, l);
          } else {
            s(k, l) = mat(l, k);
          }
        }
      }
    }
  }
  
  return S;
}
