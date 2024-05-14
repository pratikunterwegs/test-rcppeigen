
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>

// [[Rcpp::depends(RcppEigen)]]

// from https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix
// Converts an Eigen::Matrix (or expression) to Eigen::Tensor
// with dimensions specified in std::array

// template<typename Derived, typename T, auto rank>
// Eigen::Tensor<typename Derived::Scalar, rank>
// TensorCast(const Eigen::EigenBase<Derived> &matrix, const std::array<T, rank> &dims) {
//     return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>
//                 (matrix.derived().eval().data(), dims);
// }

// template<typename Derived, typename... Dims>
// auto TensorCast(const Eigen::EigenBase<Derived> &matrix, const Dims... dims) {
//     // static_assert(sizeof...(Dims) > 0, "TensorCast: sizeof... (Dims) must be larger than 0");
//     return TensorCast(matrix, std::array<Eigen::Index, sizeof...(Dims)>{dims...});
// }

auto matrix_to_tensor (const Rcpp::NumericMatrix &mat, const int& n_epi_compartments) {
    const int nrow = mat.rows();
    const int strata = mat.size() / (nrow * n_epi_compartments);

    // convert matrix to vector
    Rcpp::NumericVector mat_values(mat);

    Eigen::Tensor<double 3> tensor3d (nrow, n_epi_compartments, strata);
    // assign values manually
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < n_epi_compartments; ++j) {
            for (int k = 0; k < strata; ++k) {
                t1(i, j, k) = mat_values(i * n_epi_compartments * strata + j * strata + k);
            }
        }
    }
}

//' @export
// [[Rcpp::export()]]
void tensor_op (const Eigen::MatrixXd &mat) {
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
