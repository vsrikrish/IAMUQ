#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::cube mat_pow(int N, arma::mat A) {
  arma::cube Apow(A.n_rows, A.n_cols, N); // initialize storage for powers 0 through N-1
  // initialize first slice with I (0th power)
  Apow.slice(0) = arma::eye<arma::mat>(size(A));
  // fill each subsequent slice by A times the previous slice
  for (int s=1; s < N; ++s) {
    Apow.slice(s) = Apow.slice(s-1) * A;
  }
  return Apow; // return cube of powers
}

arma::mat sym_bind_mat(arma::cube M, arma::umat H) {
  // get matrix dimensions
  int nr_m = M.n_rows;
  int nc_m = M.n_cols;
  int nr_h = H.n_rows;
  int nc_h = H.n_cols;
  
  arma::mat S(nr_m * nr_h, nc_m * nc_h); // initialize bound matrix storage;
  // set each matrix to the appropriate cube slice
  for (int i=0; i < nr_h; ++i) {
    for (int j=0; j < nc_h; ++j) {
      int lind = H(i, j);
      if (i <= j) {
        S.submat(nr_m * i, nc_m * j, nr_m * (i+1) - 1, nc_m * (j+1) - 1) = M.slice(lind);
      } else {
        // the upper triangular matrix consists of transposes
        S.submat(nr_m * i, nc_m * j, nr_m * (i+1) - 1, nc_m * (j+1) - 1) = trans(M.slice(lind));
      }
    }
  }
  
  return S;
}

// [[Rcpp::export]]
arma::mat cov_mat(arma::mat A, arma::mat B, arma::mat D, arma::umat H) {
  // construct cube of powers of A
  int N = H.n_rows; // get number of powers to compute
  arma::cube M = mat_pow(N, A); // compute cube of powers
  M.each_slice() *= B; // multiply each slice of M by B
  arma::mat S = sym_bind_mat(M, H); // bind the appropriate submatrices together
  S += kron(arma::eye<arma::mat>(N, N), D); // add repeated elements of D along diagonal
  return S;
}

