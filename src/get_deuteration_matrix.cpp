#include <Rcpp.h>
#include <random>
#include <iterator>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

NumericMatrix get_constant_matrix(const double constant, const int nrow, const int ncol) {
    NumericVector constant_vector(nrow * ncol, constant);
    NumericMatrix constant_matrix(nrow, ncol, constant_vector.begin());
    return(constant_matrix);
}

void compare_matrix_to_scalar(const NumericMatrix& matrix, const double& value,
                              LogicalMatrix& result) {
    int N = matrix.nrow(), M = matrix.ncol();
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < M; ++j ) {
            result(i, j) = matrix(i, j) == value;
        }
    }
}

void overwrite_matrix(NumericMatrix& matrix,
                      const NumericMatrix& uniform_random,
                      const LogicalMatrix& where,
                      const NumericVector& probabilities,
                      int base) {
    int M = matrix.nrow(), N = matrix.ncol();
    for ( int i = 0; i < M; ++i ) {
        for ( int j = 0; j < N; ++j ) {
            if (where(i, j) && (uniform_random(base + i, j) <= probabilities[j])) {
                matrix(i, j) = 1.0 - matrix(i, j);
            }
        }
    }
}


//' Calculates deuteration for given timepoint
//'
//' @param initial_matrix A matrix
//' @param time_sequence vector of exchange times
//' @param hd_probs probabilities of transition HD
//' @param dh_probs probabilities of transition DH
//' @return a matrix denoting hydrogen-deuterium exchange for given timepoint.
//' @export
// [[Rcpp::export]]
NumericMatrix get_deuteration_single_timepoint(NumericMatrix initial_matrix,
                                               NumericVector time_sequence,
                                               NumericVector hd_probs,
                                               NumericVector dh_probs) {
    int M = initial_matrix.nrow();
    int N = initial_matrix.ncol();

    NumericVector v(M * N);
    v = runif(M * N);
    v.attr("dim") = Dimension(M, N);
    NumericMatrix uniform_random(Dimension(M, N));
    uniform_random = as<NumericMatrix>(v);

    LogicalMatrix HD_zeros(Dimension(M, N));
    LogicalMatrix HD_ones(Dimension(M, N));

    NumericVector result;
    for(int i = 0; i < time_sequence.length(); ++i) {
        v = runif(M *N);
        v.attr("dim") = Dimension(M, N);
        uniform_random = as<NumericMatrix>(v);

        compare_matrix_to_scalar(initial_matrix, 0.0, HD_zeros);
        compare_matrix_to_scalar(initial_matrix, 1.0, HD_ones);
        overwrite_matrix(initial_matrix, uniform_random, HD_zeros, hd_probs, 0);
        overwrite_matrix(initial_matrix, uniform_random, HD_ones, dh_probs, 0);
    }

    return(initial_matrix);
}
