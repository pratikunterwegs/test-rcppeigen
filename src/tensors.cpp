
// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

// clang-format on

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

// from
// https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix
// Converts an Eigen::Matrix (or expression) to Eigen::Tensor
// with dimensions specified in std::array

template <typename Derived, typename T, auto rank>
Eigen::Tensor<typename Derived::Scalar, rank> TensorCast(
    const Eigen::EigenBase<Derived> &matrix, const std::array<T, rank> &dims) {
  return Eigen::TensorMap<
      const Eigen::Tensor<const typename Derived::Scalar, rank>>(
      matrix.derived().eval().data(), dims);
}

template <typename Derived, typename... Dims>
auto TensorCast(const Eigen::EigenBase<Derived> &matrix, const Dims... dims) {
  // static_assert(sizeof...(Dims) > 0, "TensorCast: sizeof... (Dims) must be
  // larger than 0");
  return TensorCast(matrix, std::array<Eigen::Index, sizeof...(Dims)>{dims...});
}

auto matrix_to_tensor(const Rcpp::NumericMatrix &mat,
                      const int &n_epi_compartments) {
  const int nrow = mat.rows();
  const int strata = mat.size() / (nrow * n_epi_compartments);

  // convert matrix to vector
  Rcpp::NumericVector mat_values(mat);

  Eigen::Tensor<double, 3> tensor3d(nrow, n_epi_compartments, strata);
  // assign values manually
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < n_epi_compartments; ++j) {
      for (int k = 0; k < strata; ++k) {
        tensor3d(i, j, k) =
            mat_values(i * n_epi_compartments * strata + j * strata + k);
      }
    }
  }
}

//' @export
// [[Rcpp::export()]]
void tensor_op(const Eigen::MatrixXd &mat) {
  // const int nrow = mat.rows();
  // const int ncol = 5L;
  // const int strata = 3L;

  // create empty tensor
  Eigen::Tensor<double, 3> t1 = TensorCast(mat, 2, 2, 3);

  // // assign chip values manually
  // for (int i = 0; i < nrow; ++i) {
  //     for (int j = 0; j < ncol; ++j) {
  //         for (int k = 0; k < strata; ++k) {
  //             t1(i, j, k) = mat_values(i * ncol * strata + j * strata + k);
  //         }
  //     }
  // }

  // assign chip value
  // t3.chip(1, 1) = t3.chip(1, 1) + 0.5f;

  std::cout << t1 << "\n";
}

//' @export
// [[Rcpp::export()]]
void tensor_op2(const Rcpp::NumericVector &vec) {
  std::vector<double> vec2 = Rcpp::as<std::vector<double>>(vec);

  auto tensor = Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::ColMajor>>(
      &vec2[0], 2, 2, 3);

  std::cout << tensor.chip(0, 2) << " is the first channel \n";

  vec2[0] = 99.0;

  std::cout << tensor.chip(0, 2) << " is the first channel \n";

  tensor(0, 0, 0) = 199.0;

  std::cout << vec2[0] << "\n";
}

//' @export
// [[Rcpp::export()]]
void tensorop3() {
  /// MATRIX MULTIPLICATION OF TWO 2D TENSORS
  // Create 2 matrices using tensors of rank 2
  Eigen::Tensor<int, 2> a(2, 3);
  a.setValues({{1, 2, 3}, {4, 5, 6}});
  Eigen::Tensor<int, 2> b(3, 2);
  b.setValues({{1, 2}, {3, 4}, {5, 6}});

  // Compute the traditional matrix product
  Eigen::array<Eigen::IndexPair<int>, 1> product_dims = {
      Eigen::IndexPair<int>(1, 0)};
  Eigen::Tensor<int, 2> AB = a.contract(b, product_dims);

  // Print the resulting tensor
  std::cout << a << "\n";
  std::cout << b << "\n";
  std::cout << "Resulting tensor:\n" << AB << std::endl;

  /// ROWSUMS OF A TENSOR ALONG COLUMNS AND STRATA
  Eigen::Tensor<int, 3> x(2, 3, 3);
  x.setValues({{{0, 1, 2}, {7, 6, 5}, {8, 9, 10}},
               {{12, 13, 14}, {19, 18, 17}, {20, 21, 22}}});
  std::cout << "Resulting tensor:\n" << x << std::endl;

  Eigen::array<Eigen::Index, 3> offsets = {0, 0, 0};
  Eigen::array<Eigen::Index, 3> extents = {2, 2, 3};

  // Reduce it along the second dimension (1)...
  std::cout << "rowsumns of tensor:\n"
            << x.slice(offsets, extents).sum(Eigen::array<int, 2>({1, 2}))
            << std::endl;
}

typedef std::vector<double> state_type;

struct observer {
  std::vector<state_type> &m_states;
  std::vector<double> &m_times;

  observer(std::vector<state_type> &states,  // NOLINT
           std::vector<double> &times)       // NOLINT
      : m_states(states), m_times(times) {}

  void operator()(const state_type &x, double t) {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

//' @export
// [[Rcpp::export()]]
Rcpp::List tensor_epidemic(const double &tmax) {
  // some initial state
  state_type state{1.0, 2.0, 3.0, 4.0,  5.0,  6.0,
                   7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

  // a matrix
  Eigen::Tensor<double, 2> mat(2, 2);
  mat.setValues({{1, 2}, {3, 4}});

  // dimensions for the traditional matrix prod
  Eigen::array<Eigen::IndexPair<int>, 1> product_dims = {
      Eigen::IndexPair<int>(1, 0)};

  // an odeint-compatible functionobject
  struct model {
    const double beta = 1.3 / 7.0;
    const double gamma = 1.0 / 7.0;
    Eigen::Tensor<double, 2> matrix;

    explicit model(Eigen::Tensor<double, 2> matrix) : matrix(matrix) {}

    void operator()(state_type &x, state_type &dx, const double &t) {
      // map a tensor to the vector
      auto x_tensor =
          Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::ColMajor>>(&x[0], 2,
                                                                      2, 3);
      auto dx_tensor =
          Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::ColMajor>>(&dx[0], 2,
                                                                      2, 3);

      // implement matrix mult on slices
      Eigen::Tensor<double, 2> a(2, 2);
      a.setValues({{9.0, 1.0}, {9.0, 1.0}});

      // Compute the traditional matrix product
      Eigen::array<Eigen::IndexPair<int>, 1> product_dims = {
          Eigen::IndexPair<int>(1, 0)};

      for (size_t i = 0; i < 3; i++) {
        dx_tensor.chip(i, 2).chip(0, 1) =
            a.contract(dx_tensor.chip(i, 2).chip(0, 1) + 1.0, product_dims);
      }

      dx_tensor(0, 0, 0) = dx_tensor(0, 0, 0) + (99999.0 * t);
    }
  };

  // define model
  model this_model(mat);

  // prepare storage containers for the observer
  std::vector<state_type> x_vec;  // is a vector of vectors
  std::vector<double> times;

  // a controlled stepper for constant step sizes
  boost::numeric::odeint::runge_kutta4<state_type> stepper;

  // run the function without assignment
  boost::numeric::odeint::integrate_const(stepper, this_model, state, 0.0, tmax,
                                          0.1, observer(x_vec, times));

  return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x_vec),
                            Rcpp::Named("time") = Rcpp::wrap(times));
}