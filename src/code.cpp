#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//' Calculate squares.
//'
//' @param x A vector of numbers to square.
//' @return Square of x.
//'
// [[Rcpp::export]]
Eigen::VectorXd square_rcppeigen_internal(Eigen::VectorXd x) {
  Eigen::VectorXd y = x.array().square();
  return(y);
}
