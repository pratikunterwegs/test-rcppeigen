
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>

// [[Rcpp::depends(RcppEigen)]]

// from https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix
// Converts an Eigen::Matrix (or expression) to Eigen::Tensor
// with dimensions specified in std::array

template<typename Derived, typename T, auto rank>
Eigen::Tensor<typename Derived::Scalar, rank>
TensorCast(const Eigen::EigenBase<Derived> &matrix, const std::array<T, rank> &dims) {
    return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>
                (matrix.derived().eval().data(), dims);
}

template<typename Derived, typename... Dims>
auto TensorCast(const Eigen::EigenBase<Derived> &matrix, const Dims... dims) {
    // static_assert(sizeof...(Dims) > 0, "TensorCast: sizeof... (Dims) must be larger than 0");
    return TensorCast(matrix, std::array<Eigen::Index, sizeof...(Dims)>{dims...});
}

auto matrix_to_tensor (const Rcpp::NumericMatrix &mat, const int& n_epi_compartments) {
    const int nrow = mat.rows();
    const int strata = mat.size() / (nrow * n_epi_compartments);

    // convert matrix to vector
    Rcpp::NumericVector mat_values(mat);

    Eigen::Tensor<double, 3> tensor3d (nrow, n_epi_compartments, strata);
    // assign values manually
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < n_epi_compartments; ++j) {
            for (int k = 0; k < strata; ++k) {
                tensor3d(i, j, k) = mat_values(i * n_epi_compartments * strata + j * strata + k);
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

//' @export
// [[Rcpp::export()]]
void tensor_op2 (const Rcpp::NumericVector &vec) {
    std::vector<double> vec2 = Rcpp::as<std::vector<double> > (vec);

    auto tensor = Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::ColMajor>> (&vec2[0], 2, 2, 3);

    std::cout << tensor.chip(0, 2) << " is the first channel \n";

    vec2[0] = 99.0;

    std::cout << tensor.chip(0, 2) << " is the first channel \n";

    std::cout << tensor.chip(0, 1) * tensor.chip(0, 2)  << "\nare the first cols \n";
}

//' @export
// [[Rcpp::export()]]
void tensorop3 () {

    /// MATRIX MULTIPLICATION OF TWO 2D TENSORS
    // Create 2 matrices using tensors of rank 2
    Eigen::Tensor<int, 2> a(2, 3);
    a.setValues({{1, 2, 3}, {4, 5, 6}});
    Eigen::Tensor<int, 2> b(3, 2);
    b.setValues({{1, 2}, {3, 4}, {5, 6}});

    // Compute the traditional matrix product
    Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>(1, 0) };
    Eigen::Tensor<int, 2> AB = a.contract(b, product_dims);

    // Print the resulting tensor
    std::cout << a << "\n";
    std::cout << b << "\n";
    std::cout << "Resulting tensor:\n" << AB << std::endl;

    /// ROWSUMS OF A TENSOR ALONG COLUMNS AND STRATA
    Eigen::Tensor<int, 3> x (2, 3, 3);
    x.setValues({{{0, 1, 2},
              {7, 6, 5},
              {8, 9, 10}},
             {{12, 13, 14},
              {19, 18, 17},
              {20, 21, 22}}});
    std::cout << "Resulting tensor:\n" << x << std::endl;

    Eigen::array<Eigen::Index, 3> offsets = {0, 0, 0};
    Eigen::array<Eigen::Index, 3> extents = {2, 2, 3};

    // Reduce it along the second dimension (1)...
    std::cout << "rowsumns of tensor:\n" << x.slice(offsets, extents).sum(Eigen::array<int, 2>({1, 2})) << std::endl;
}
