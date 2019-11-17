#include "points_index.h"

// from dodgr/src/sf_as_network.cpp

struct OnePointIndex : public RcppParallel::Worker
{
    const Rcpp::NumericVector xy_x, xy_y, pt_x, pt_y;
    const size_t nxy;
    RcppParallel::RVector <int> index;

    // constructor
    OnePointIndex (
            const Rcpp::NumericVector xy_x_in,
            const Rcpp::NumericVector xy_y_in,
            const Rcpp::NumericVector pt_x_in,
            const Rcpp::NumericVector pt_y_in,
            const size_t nxy_in,
            Rcpp::IntegerVector index_in) :
        xy_x (xy_x_in), xy_y (xy_y_in), pt_x (pt_x_in), pt_y (pt_y_in),
        nxy (nxy_in), index (index_in)
    {
    }

    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        for (std::size_t i = begin; i < end; i++)
        {
            long int li = static_cast <long int> (i);
            double dmin = INFINITE_DOUBLE;
            int jmin = INFINITE_INT;
            for (int j = 0; j < static_cast <int> (nxy); j++)
            {
                double dij = (xy_x [j] - pt_x [li]) * (xy_x [j] - pt_x [li]) +
                    (xy_y [j] - pt_y [li]) * (xy_y [j] - pt_y [li]);
                if (dij < dmin)
                {
                    dmin = dij;
                    jmin = j;
                }
            }
            index [i] = jmin;
        }
    }
                                   
};

//' rcpp_points_index_par
//'
//' Get index of nearest vertices to list of points
//'
//' @param graph Rcpp::DataFrame containing the graph
//' @param pts Rcpp::DataFrame containing the routing points
//'
//' @return 0-indexed Rcpp::NumericVector index into graph of nearest points
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_points_index_par (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts)
{
    Rcpp::NumericVector ptx = pts ["x"];
    Rcpp::NumericVector pty = pts ["y"];

    Rcpp::NumericVector vtx = xy ["x"];
    Rcpp::NumericVector vty = xy ["y"];

    size_t npts = static_cast <size_t> (pts.nrow ()),
           nxy = static_cast <size_t> (xy.nrow ());

    //Rcpp::IntegerVector index (n, Rcpp::IntegerVector::get_na ());
    Rcpp::IntegerVector index (npts);
    // Create parallel worker
    OnePointIndex one_pt_indx (vtx, vty, ptx, pty, nxy, index);

    RcppParallel::parallelFor (0, npts, one_pt_indx);

    return index;
}
