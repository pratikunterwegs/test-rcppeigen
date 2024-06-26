#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>

#include <./unsupported/Eigen/CXX11/Tensor>

// [[Rcpp::depends(RcppEigen)]]
using namespace std;
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Rcpp::List epi_spread(
    const Eigen::MatrixXd &contact_matrix,
    const Eigen::VectorXd &demography_vector,
    const Eigen::MatrixXd &p_susceptibility,  // risk groups
    const Eigen::MatrixXd &susceptibility     // susc of risk grp?
) {
  // count number of risk groups
  const int n_susc_groups = p_susceptibility.cols();

  // make single column matrix from prop_suscep data,
  // prop_suscep is the prob(suscep) per demography group
  Eigen::MatrixXd psusc = p_susceptibility;
  const Eigen::VectorXd lps(
      Eigen::Map<Eigen::VectorXd>(psusc.data(), psusc.size()));

  // replicate demography vector by the number of risk groups
  // and multiply it by the prop_suscep values
  const Eigen::VectorXd demography_vector_ =
      demography_vector.replicate(n_susc_groups, 1).array() * lps.array();

  const Eigen::MatrixXd contact_matrix_ =
      contact_matrix.replicate(n_susc_groups, n_susc_groups);

  // unroll the risk level matrix
  Eigen::MatrixXd susc = susceptibility;
  const Eigen::VectorXd rm(
      Eigen::Map<Eigen::VectorXd>(susc.data(), susc.size()));

  return Rcpp::List::create(
      Rcpp::Named("contact_matrix") = contact_matrix_,
      Rcpp::Named("demography_vector") = demography_vector_,
      Rcpp::Named("susceptibility") = rm);
}

template <typename T>
Rcpp::NumericVector copyFromTensor(const T &x) {
  int n = x.size();
  NumericVector ans(n);
  IntegerVector dim(x.NumDimensions);
  for (int i = 0; i < x.NumDimensions; ++i) {
    dim[i] = x.dimension(i);
  }
  memcpy((double *)ans.begin(), (double *)x.data(), n * sizeof(double));
  ans.attr("dim") = dim;
  return ans;
}

// [[Rcpp::export]]
NumericVector getTensor() {
  Eigen::Tensor<double, 2> x(2, 3);
  x.setValues({{1, 2, 3}, {6, 5, 4}});
  Eigen::array<int, 1> dims({1 /* dimension to reduce */});
  return NumericVector(x.size());
}
