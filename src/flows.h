#pragma once

#include <memory>
#include <vector>
#include <algorithm> // std::fill, std::reverse
#include <iostream>
#include <fstream>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include "pathfinders.h"

class DGraph;
class PathFinder;

//----------------------------
//----- functions in flows.cpp
//----------------------------

size_t get_chunk_size (const size_t nfrom);

Rcpp::NumericVector rcpp_flows_aggregate_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        Rcpp::NumericMatrix flows,
        const bool norm_sums,
        const double tol);

Rcpp::NumericVector rcpp_flows_disperse_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::NumericVector k,
        Rcpp::NumericVector flows,
        const double &tol);

Rcpp::NumericVector rcpp_flows_si (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        Rcpp::NumericVector kvec,
        Rcpp::NumericVector dens_from,
        Rcpp::NumericVector dens_to,
        const bool norm_sums,
        const double tol);

Rcpp::NumericVector rcpp_flows_exposure (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::NumericVector k,
        Rcpp::NumericVector dens,
        const double &tol);

