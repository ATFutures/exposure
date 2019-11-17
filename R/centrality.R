#' exposure_centrality
#'
#' Calculate time-based betweenness centrality for a street network, as routed
#' either for vehicular or pedestrian transport.
#'
#' @param graph 'data.frame' or equivalent object representing the network
#' graph (see Details)
#' @param pedestrian If `TRUE`, calculate pedestrian-weighted centrality,
#' otherwise centrality for vehicular (motorcar) transport.
#' @return Modified version of graph with additonal 'centrality' column added.
#'
#' @export
exposure_centrality <- function (graph, pedestrian = TRUE)
{
    if ("centrality" %in% names (graph))
        warning ("graph already has a 'centrality' column; ",
                  "this will be overwritten")

    dist_threshold <- .Machine$double.xmax
    edges <- TRUE # hard-coded for edge-based centrality
    heap <- "BHeap"

    gr_cols <- dodgr_graph_cols (graph)
    vert_map <- make_vert_map (graph, gr_cols)
    graph2 <- convert_graph (graph, gr_cols)

    # final '0' is for sampling calculation to estimate speed - non-zero values
    # used only in 'estimate_centrality_time'
    centrality <- rcpp_centrality (graph2, vert_map, heap, dist_threshold, edges, 0)

    graph$centrality <- centrality

    return (graph)
}
