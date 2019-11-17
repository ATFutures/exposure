#' exposure_centrality
#'
#' Calculate time-based betweenness centrality for a street network, as routed
#' either for vehicular or pedestrian transport.
#'
#' @param net An `sc`-class street network extracted with
#' `dodgr::dodgr_streetnet_sc`.
#' @param dispersal Range in metres over which vehicular emissions are
#' effectively dispersed.
#' @return A `dodgr` network with additional columns of 'centrality-veh' for
#' vehicular centrality; 'centrality-ped' for pedestrian centrality; and
#' 'exposure' for relative risk from exposure of pedestrians to vehicular
#' emissions.
#'
#' @export
exposure_centrality <- function (net, dispersal = 20)
{
    net_p <- centrality1 (net, wt_profile = "foot")
    net_v <- centrality1 (net, wt_profile = "motorcar")
    disperse_emissions (net_p, net_v, dispersal = dispersal)
}

centrality1 <- function (net, wt_profile = "foot")
{
    message (cli::col_blue (cli::symbol$pointer),
             cli::col_green ("  Preparing streetnet for ",
                             wt_profile, " weighting"), appendLF = FALSE)

    dodgr::dodgr_cache_off ()
    net <- dodgr::weight_streetnet (net, wt_profile = wt_profile)
    net$d <- net$time
    net$d_weighted <- net$time_weighted
    net_c <- dodgr::dodgr_contract_graph (net)

    dist_threshold <- .Machine$double.xmax
    edges <- TRUE # hard-coded for edge-based centrality

    gr_cols <- dodgr_graph_cols (net_c)
    vert_map <- make_vert_map (net_c, gr_cols)
    graph <- convert_graph (net_c, gr_cols)

    message ("\r", cli::col_green (cli::symbol$tick,
             "  Prepared streetnet for ", wt_profile, " weighting "))
    message (cli::col_blue (cli::symbol$pointer),
             cli::col_green ("  Calculating centrality"), appendLF = FALSE)

    # final '0' is for sampling calculation to estimate speed - non-zero values
    # used only in 'estimate_centrality_time'
    centrality <- rcpp_centrality (graph, vert_map, dist_threshold, edges, 0)

    message ("\r", cli::col_green (cli::symbol$tick,
             "  Calculated centrality "))

    net_c$centrality <- centrality
    dodgr::dodgr_uncontract_graph (net_c)
}

# Disperse relative emissions from net_v back on to net_p, noting that these two
# networks will generally differ due to different weighting profiles
disperse_emissions <- function (net_p, net_v, dispersal = 20)
{
    message (cli::col_blue (cli::symbol$pointer),
             cli::col_green ("  Dispersing emissions between networks"),
             appendLF = FALSE)

    v <- rbind (dodgr::dodgr_vertices (net_v),
                dodgr::dodgr_vertices (net_p))
    v <- v [match (unique (v$id), v$id), ]
    # append net_v centrality values to v:
    indxf <- match (v$id, net_v$.vx0)
    indxt <- match (v$id, net_v$.vx1)
    indxf [is.na (indxf)] <- indxt [is.na (indxf)]
    indxt [is.na (indxt)] <- indxf [is.na (indxt)]
    v$centrality <- apply (cbind (net_v$centrality [indxf],
                                  net_v$centrality [indxt]), 1, max)

    # match any NA values from IDs in net_p but not in net_v to closest net_v
    # vertices
    index <- which (is.na (v$centrality))
    vv <- dodgr::dodgr_vertices (net_v)
    ids <- vv$id [dodgr::match_points_to_graph (vv, v [index, c ("x", "y")])]
    v$centrality [index] <- v$centrality [match (ids, v$id)]

    # disperse centrality values via Gaussian kernel
    # warnings are issued here if any (x, y) are duplicated
    xy <- suppressWarnings (spatstat::ppp (v$x, v$y, range (v$x), range (v$y)))
    dist_to_lonlat_range <- function (verts, d = 20)
    {
        xy0 <- c (mean (verts$x), mean (verts$y))
        names (xy0) <- c ("x", "y")
        minf <- function (a, xy0) { abs (geodist::geodist (xy0, xy0 + a) - d) }
        stats::optimise (minf, c (0, 0.1), xy0)$minimum
    }
    sig <- dist_to_lonlat_range (v, d = dispersal)
    d <- spatstat::density.ppp (xy, weights = v$centrality, sigma = sig,
                                at = "points")
    vsum <- sum (v$centrality)
    d <- d * vsum / sum (d)
    indx <- which (d > v$centrality)
    v$centrality [indx] <- d [indx]
    v$centrality <- v$centrality * vsum / sum (v$centrality)

    # Then map those vertex values back onto net_p
    indxf <- match (net_p$.vx0, v$id)
    indxt <- match (net_p$.vx1, v$id)
    fmax <- apply (cbind (v$centrality [indxf], v$centrality [indxt]), 1, max)
    fsum <- sum (net_v$centrality)
    fmax <- fmax * fsum / sum (net_p$centrality)
    indx <- which (fmax > net_p$centrality)
    # That again imbalances the centrality of the graph, so needs to be standardised
    # back to original sum again
    fsum <- sum (net_p$centrality)
    net_p$exposure <- 0
    net_p$exposure [indx] <- fmax [indx]
    net_p$exposure <- net_p$exposure * fsum / sum (net_p$exposure)

    message ("\r", cli::col_green (cli::symbol$tick,
             "  Dispersed emissions between networks "))

    return (net_p)
}
