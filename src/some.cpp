
#include <Rcpp.h>
#include <RcppEigen.h>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <fstream>
#include <iostream>
#include <utility>

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

//[ stiff_system_definition
typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

struct stiff_system {
  void operator()(const vector_type &x, vector_type &dxdt, double /* t */) {
    dxdt[0] = -101.0 * x[0] - 100.0 * x[1];
    dxdt[1] = x[0];
  }
};

struct stiff_system_jacobi {
  void operator()(const vector_type & /* x */, matrix_type &J,
                  const double & /* t */, vector_type &dfdt) {
    J(0, 0) = -101.0;
    J(0, 1) = -100.0;
    J(1, 0) = 1.0;
    J(1, 1) = 0.0;
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
  }
};
//]

//' @export
// [[Rcpp::export]]
Rcpp::List somefun() {
  vector_type x(2, 1.0);
  size_t num_of_steps = integrate_const(
      make_dense_output<rosenbrock4<double> >(1.0e-6, 1.0e-6),
      make_pair(stiff_system(), stiff_system_jacobi()), x, 0.0, 50.0, 0.01,
      Rcpp::Rcout << phoenix::arg_names::arg2 << " "
                  << phoenix::arg_names::arg1[0] << "\n");
  //]

  // using rk_dopri5
  vector_type x2(2, 1.0);
  size_t num_of_steps2 = integrate_const(
      make_dense_output<runge_kutta_dopri5<vector_type> >(1.0e-6, 1.0e-6),
      stiff_system(), x2, 0.0, 50.0, 0.01,
      Rcpp::Rcout << phoenix::arg_names::arg2 << " "
                  << phoenix::arg_names::arg1[0] << "\n");
  //]

  return Rcpp::List::create(num_of_steps, num_of_steps2);
}
