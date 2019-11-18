#' exposure
#'
#' Calculate time-based betweenness centrality for a street network, as routed
#' either for vehicular or pedestrian transport.
#'
#' @param net_p A \pkg{dodgr} network weighted for pedestrian transport
#' @param net_v A \pkg{dodgr} network weighted for vehicular transport with an
#' additional `centrality` column appended by \link{exposure_centrality}.
#' @param dispersal Range in metres over which vehicular emissions are
#' effectively dispersed.
#' @return A `dodgr` network with additional columns of 'centrality-veh' for
#' vehicular centrality; 'centrality-ped' for pedestrian centrality; and
#' 'exposure' for relative risk from exposure of pedestrians to vehicular
#' emissions.
#'
#' @export
exposure <- function (net_p, net_v, dispersal = 20)
{
    if ("centrality" %in% names (net_p) &&
        !"centrality" %in% names (net_v))
        stop ("parameters are 'net_p' then 'net_v'")

    disperse_emissions (net_p, net_v, dispersal = dispersal)
}

#' exposure_centrality
#'
#' Calculate centrality for a street network
#'
#' @param net A \pkg{dodgr} network weighted for vehicular transport (that is,
#' with `wt_profile = "motorcar"`).
#' @return Modified version of same net with additional "centrality" column
#' appended.
#' @export
exposure_centrality <- function (net)
{
    message (cli::col_blue (cli::symbol$pointer),
             cli::col_green ("  Preparing network "), appendLF = FALSE)

    dodgr::dodgr_cache_off ()
    net$d <- net$time
    net$d_weighted <- net$time_weighted
    net_c <- dodgr::dodgr_contract_graph (net)

    message ("\r", cli::col_green (cli::symbol$tick,
             "  Prepared network "))
    message (cli::col_blue (cli::symbol$pointer),
             cli::col_green ("  Calculating centrality"), appendLF = FALSE)

    # final '0' is for sampling calculation to estimate speed - non-zero values
    # used only in 'estimate_centrality_time'
    net_c <- dodgr::dodgr_centrality  (net_c, contract = FALSE,
                                       edges = TRUE)

    message ("\r", cli::col_green (cli::symbol$tick,
             "  Calculated centrality "))

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
    if ("from_id" %in% names (net_v))
    {
        indxf <- match (v$id, net_v$from_id)
        indxt <- match (v$id, net_v$to_id)
    } else {
        indxf <- match (v$id, net_v$.vx0)
        indxt <- match (v$id, net_v$.vx1)
    }
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
    v$exposure <- spatstat::density.ppp (xy, weights = v$centrality,
                                         sigma = sig, at = "points")

    # Then map those vertex values back onto net_p
    if ("from_id" %in% names (net_v))
    {
        indxf <- match (net_p$from_id, v$id)
        indxt <- match (net_p$to_id, v$id)
    } else {
        indxf <- match (net_p$.vx0, v$id)
        indxt <- match (net_p$.vx1, v$id)
    }
    fmax <- apply (cbind (v$exposure [indxf], v$exposure [indxt]), 1, max)
    net_p$exposure <- fmax

    message ("\r", cli::col_green (cli::symbol$tick,
             "  Dispersed emissions between networks "))

    return (net_p)
}
