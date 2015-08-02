// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// quantileCpp
NumericVector quantileCpp(NumericVector x, NumericVector probs);
RcppExport SEXP D3M_quantileCpp(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    __result = Rcpp::wrap(quantileCpp(x, probs));
    return __result;
END_RCPP
}
// wasserCpp
double wasserCpp(NumericVector x, NumericVector y, int paranum, int q);
RcppExport SEXP D3M_wasserCpp(SEXP xSEXP, SEXP ySEXP, SEXP paranumSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type paranum(paranumSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    __result = Rcpp::wrap(wasserCpp(x, y, paranum, q));
    return __result;
END_RCPP
}
// wasserCpp_mat
NumericVector wasserCpp_mat(NumericMatrix xMat, NumericMatrix yMat, int paranum, int q);
RcppExport SEXP D3M_wasserCpp_mat(SEXP xMatSEXP, SEXP yMatSEXP, SEXP paranumSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type xMat(xMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< int >::type paranum(paranumSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    __result = Rcpp::wrap(wasserCpp_mat(xMat, yMat, paranum, q));
    return __result;
END_RCPP
}
// randomShuffle
Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a);
RcppExport SEXP D3M_randomShuffle(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type a(aSEXP);
    __result = Rcpp::wrap(randomShuffle(a));
    return __result;
END_RCPP
}
// permCpp
List permCpp(NumericMatrix casesMat, NumericMatrix controlMat, NumericVector d, int bsn, int qn, int q);
RcppExport SEXP D3M_permCpp(SEXP casesMatSEXP, SEXP controlMatSEXP, SEXP dSEXP, SEXP bsnSEXP, SEXP qnSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type casesMat(casesMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type controlMat(controlMatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type bsn(bsnSEXP);
    Rcpp::traits::input_parameter< int >::type qn(qnSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    __result = Rcpp::wrap(permCpp(casesMat, controlMat, d, bsn, qn, q));
    return __result;
END_RCPP
}
