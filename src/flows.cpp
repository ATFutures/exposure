
#include "run_sp.h"
#include "flows.h"

#include "dgraph.h"
#include "heaps/heap_lib.h"

template <typename T>
void inst_graph (std::shared_ptr<DGraph> g, unsigned int nedges,
        const std::map <std::string, unsigned int>& vert_map,
        const std::vector <std::string>& from,
        const std::vector <std::string>& to,
        const std::vector <T>& dist,
        const std::vector <T>& wt)
{
    for (unsigned int i = 0; i < nedges; ++i)
    {
        unsigned int fromi = vert_map.at(from [i]);
        unsigned int toi = vert_map.at(to [i]);
        g->addNewEdge (fromi, toi, dist [i], wt [i], i);
    }
}


struct OneExposure : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const Rcpp::NumericVector dens;
    const std::vector <std::string> vert_name;
    const std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    size_t nverts; // can't be const because of reinterpret cast
    size_t nedges;
    const double k;
    const double tol;
    std::shared_ptr <DGraph> g;

    std::vector <double> output;

    // Constructor 1: The main constructor
    OneExposure (
            const Rcpp::IntegerVector fromi,
            const Rcpp::NumericVector dens_in,
            const std::vector <std::string>  vert_name_in,
            const std::unordered_map <std::string, unsigned int> verts_to_edge_map_in,
            const size_t nverts_in,
            const size_t nedges_in,
            const double k_in,
            const double tol_in,
            const std::shared_ptr <DGraph> g_in) :
        dp_fromi (fromi), dens (dens_in), vert_name (vert_name_in),
        verts_to_edge_map (verts_to_edge_map_in),
        nverts (nverts_in), nedges (nedges_in), k (k_in),
        tol (tol_in), g (g_in), output ()
    {
        output.resize (nverts, 0.0);
    }

    // Constructor 2: The Split constructor
    OneExposure (
            const OneExposure& oneExposure,
            RcppParallel::Split) :
        dp_fromi (oneExposure.dp_fromi), dens (oneExposure.dens),
        vert_name (oneExposure.vert_name),
        verts_to_edge_map (oneExposure.verts_to_edge_map),
        nverts (oneExposure.nverts), nedges (oneExposure.nedges),
        k (oneExposure.k), tol (oneExposure.tol),
        g (oneExposure.g), output ()
    {
        output.resize (nverts, 0.0);
    }


    // Parallel function operator
    void operator() (size_t begin, size_t end)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                    *run_sp::getHeapImpl (), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <int> prev (nverts);

        for (size_t i = begin; i < end; i++) // over the from vertices
        {
            // translate k-value to distance limit based on tol
            // exp(-d / k) = tol -> d = -k * log (tol)
            double dlim = -k * log (tol);

            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);
            std::fill (prev.begin (), prev.end (), INFINITE_INT);

            const unsigned int from_i = static_cast <unsigned int> (dp_fromi [i]);
            d [from_i] = w [from_i] = 0.0;

            pathfinder->DijkstraLimit (d, w, prev, from_i, dlim);

            // get sums tracing back from terminal vertices
            std::vector <bool> has_prev (nverts, false);
            std::unordered_set <int> prev_set;
            for (unsigned int j = 0; j < nverts; j++)
            {
                if (w [j] < dlim)
                    prev [i] = -1;
                else
                {
                    has_prev [i] = true;
                    if (prev [i] > -1)
                        prev_set.emplace (prev [i]);
                }
            }
            double expsum = 0.0, n = 0.0;
            for (unsigned int j = 0; j < nverts; j++)
            {
                if (has_prev [j] && prev_set.find (j) == prev_set.end ())
                {
                    expsum += d [i];
                    n += 1.0;
                }
                    //output [i] += d [i];
            }
            // calculate average exposure from that point:
            if (n > 0.0)
                output [i] += expsum / n;
        } // end for i
    } // end parallel function operator

    void join (const OneExposure &rhs)
    {
        for (size_t i = 0; i < nverts; i++)
            output [i] += rhs.output [i];
    }
};

// RcppParallel jobs can be chunked to a specified "grain size"; see
// https://rcppcore.github.io/RcppParallel/#grain_size
// This function determines chunk size such that there are at least 100 chunks
// for a given `nfrom`.
size_t get_chunk_size (const size_t nfrom)
{
    size_t chunk_size;

    if (nfrom > 1000)
        chunk_size = 100;
    else if (nfrom > 100)
        chunk_size = 10;
    else
        chunk_size = 1;

    return chunk_size;
}


//' rcpp_flows_exposure
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param fromi Index into vert_map_in of vertex numbers
//' @param k Coefficient of (current proof-of-principle-only) exponential
//' distance decay function.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_flows_exposure (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        const double k,
        Rcpp::NumericVector dens,
        const double &tol)
{
    Rcpp::NumericVector id_vec;
    const size_t nfrom = static_cast <size_t> (fromi.size ());

    std::vector <std::string> from = graph [".vx0"];
    std::vector <std::string> to = graph [".vx1"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["d_weighted"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::vector <std::string> vert_name = vert_map_in ["vert"];
    std::vector <unsigned int> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, unsigned int> vert_map_i;
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    run_sp::make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    // Create parallel worker
    OneExposure oneExposure (fromi, dens, vert_name, verts_to_edge_map,
            nverts, nedges, k, tol, g);

    size_t chunk_size = get_chunk_size (nfrom);
    RcppParallel::parallelReduce (0, nfrom, oneExposure, chunk_size);

    return Rcpp::wrap (oneExposure.output);
}


