#' exposure
#'
#' Disperse flows throughout a network based on a input vectors of origin points
#' and associated densities, and calculate resultant relative risk of exposure
#' to urban pollutants.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph, with one column quantifying pollutant concentrations throughout the
#' network.
#' @param from Vector or matrix of points **from** which aggregate dispersed
#' flows are to be calculated (see Details)
#' @param dens Vectors of densities correponsing to the `from` points
#' @param k Width coefficient of exponential diffusion function defined as
#' `exp(-d/k)`, in units of distance column of `graph` (metres by default).
#' @param tol Relative tolerance below which dispersal is considered to have
#' finished. This parameter can generally be ignored; if in doubt, its effect
#' can be removed by setting `tol = 0`.
#' @return Modified version of graph with additonal `exposure` column added.
#'
#' @export
exposure <- function (graph, from, dens, k = 500, tol = 1e-12)
{
    if (any (is.na (dens))) {
        dens [is.na (dens)] <- 0
    }
    k <- k [1]

    heap <- "BHeap"
    g <- prepare_graph (graph, from)

    if (!is.matrix (dens))
        dens <- as.matrix (dens)

    f <- rcpp_flows_disperse_par (g$graph, g$vert_map, g$from_index,
                                  k, dens, tol, heap)
    graph$exposure <- f

    return (graph)
}


# transform input graph and (from, to) arguments to standard forms for passing
# to C++ routines
prepare_graph <- function (graph, from, to)
{
    if ("flow" %in% names (graph))
        warning ("graph already has a 'flow' column; ",
                  "this will be overwritten")

    gr_cols <- dodgr_graph_cols (graph)
    vert_map <- make_vert_map (graph, gr_cols)

    # change from and to just to check conformity
    tp <- attr (graph, "turn_penalty")
    tp <- ifelse (is.null (tp), 0, tp)

    # remove any routing points not in edge start nodes:
    from <- nodes_arg_to_pts (from, graph)
    if (methods::is (graph, "dodgr_streetnet_sc") & tp > 0)
        from <- remap_verts_with_turn_penalty (graph, from, from = TRUE)
    from <- from [from %in% graph [[gr_cols$from]] ]
    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based

    to_index <- NULL
    if (!missing (to))
    {
        # remove any routing points not in edge end nodes:
        to <- nodes_arg_to_pts (to, graph)
        if (methods::is (graph, "dodgr_streetnet_sc") & tp > 0)
            to <- remap_verts_with_turn_penalty (graph, to, from = FALSE)
        to <- to [to %in% graph [[gr_cols$to]] ]
        index_id <- get_index_id_cols (graph, gr_cols, vert_map, to)
        to_index <- index_id$index - 1 # 0-based
    }

    graph2 <- convert_graph (graph, gr_cols)

    list (graph = graph2, vert_map = vert_map,
          from_index = from_index, to_index = to_index)
}

nodes_arg_to_pts <- function (nodes, graph)
{
    if (!is.matrix (nodes))
        nodes <- as.matrix (nodes)
    if (ncol (nodes) == 2)
    {
        verts <- dodgr::dodgr_vertices (graph)
        nodes <- verts$id [dodgr::match_points_to_graph (verts, nodes)]
    }
    return (nodes)
}
