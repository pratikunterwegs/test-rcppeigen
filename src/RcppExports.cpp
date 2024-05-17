// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// epi_spread
Rcpp::List epi_spread(const Eigen::MatrixXd& contact_matrix, const Eigen::VectorXd& demography_vector, const Eigen::MatrixXd& p_susceptibility, const Eigen::MatrixXd& susceptibility);
RcppExport SEXP _testrcppeigen_epi_spread(SEXP contact_matrixSEXP, SEXP demography_vectorSEXP, SEXP p_susceptibilitySEXP, SEXP susceptibilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type contact_matrix(contact_matrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type demography_vector(demography_vectorSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type p_susceptibility(p_susceptibilitySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type susceptibility(susceptibilitySEXP);
    rcpp_result_gen = Rcpp::wrap(epi_spread(contact_matrix, demography_vector, p_susceptibility, susceptibility));
    return rcpp_result_gen;
END_RCPP
}
// getTensor
NumericVector getTensor();
RcppExport SEXP _testrcppeigen_getTensor() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getTensor());
    return rcpp_result_gen;
END_RCPP
}
// somefun
Rcpp::List somefun();
RcppExport SEXP _testrcppeigen_somefun() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(somefun());
    return rcpp_result_gen;
END_RCPP
}
// tensor_op
void tensor_op(const Eigen::MatrixXd& mat);
RcppExport SEXP _testrcppeigen_tensor_op(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type mat(matSEXP);
    tensor_op(mat);
    return R_NilValue;
END_RCPP
}
// tensor_op2
void tensor_op2(const Rcpp::NumericVector& vec);
RcppExport SEXP _testrcppeigen_tensor_op2(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    tensor_op2(vec);
    return R_NilValue;
END_RCPP
}
// tensorop3
void tensorop3();
RcppExport SEXP _testrcppeigen_tensorop3() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    tensorop3();
    return R_NilValue;
END_RCPP
}
// tensor_epidemic
Rcpp::List tensor_epidemic(const double& tmax);
RcppExport SEXP _testrcppeigen_tensor_epidemic(SEXP tmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type tmax(tmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(tensor_epidemic(tmax));
    return rcpp_result_gen;
END_RCPP
}
// aging
void aging();
RcppExport SEXP _testrcppeigen_aging() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    aging();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_testrcppeigen_epi_spread", (DL_FUNC) &_testrcppeigen_epi_spread, 4},
    {"_testrcppeigen_getTensor", (DL_FUNC) &_testrcppeigen_getTensor, 0},
    {"_testrcppeigen_somefun", (DL_FUNC) &_testrcppeigen_somefun, 0},
    {"_testrcppeigen_tensor_op", (DL_FUNC) &_testrcppeigen_tensor_op, 1},
    {"_testrcppeigen_tensor_op2", (DL_FUNC) &_testrcppeigen_tensor_op2, 1},
    {"_testrcppeigen_tensorop3", (DL_FUNC) &_testrcppeigen_tensorop3, 0},
    {"_testrcppeigen_tensor_epidemic", (DL_FUNC) &_testrcppeigen_tensor_epidemic, 1},
    {"_testrcppeigen_aging", (DL_FUNC) &_testrcppeigen_aging, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_testrcppeigen(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
