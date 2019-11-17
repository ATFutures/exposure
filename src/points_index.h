#include <string>
#include <cmath>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

constexpr float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
constexpr int INFINITE_INT =  std::numeric_limits<int>::max ();

Rcpp::IntegerVector rcpp_points_index_par (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts);
