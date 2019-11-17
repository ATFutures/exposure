#' dodgr_graph_cols
#'
#' Identify the essential columns of the graph table (data.frame, tibble,
#' whatever) to be analysed in the C++ routines.
#'
#' @param graph A `data.frame` containing the edges of the graph
#' @return A list of column numbers of `edge_id`, `from`,
#' `to`, `d`, `w`, `time`, `xfr`, `yfr`, `xto`, `yto`, and `component`, some of
#' which may be NA.
#'
#' @noRd
dodgr_graph_cols <- function (graph)
{
    nms <- names (graph)
    component <- grep ("comp", nms) %>% null_to_na ()
    if (methods::is (graph, "dodgr_streetnet") &
        !methods::is (graph, "dodgr_streetnet_sc") & ncol (graph) >= 11)
    {
        # columns are always identically structured
        edge_id <- which (nms == "edge_id") %>% null_to_na ()
        fr_col <- which (nms == "from_id") %>% null_to_na ()
        to_col <- which (nms == "to_id") %>% null_to_na ()
        d_col <- which (nms == "d")
        w_col <- which (nms == "d_weighted")

        xfr <- which (nms == "from_lon")
        if (length (xfr) == 0) xfr <- NA
        yfr <- which (nms == "from_lat")
        if (length (yfr) == 0) yfr <- NA
        xto <- which (nms == "to_lon")
        if (length (xto) == 0) xto <- NA
        yto <- which (nms == "to_lat")
        if (length (yto) == 0) yto <- NA
    } else
    {
        edge_id <- grep ("edge_id|edge_$", nms) %>% null_to_na ()

        d_col <- find_d_col (graph)
        w_col <- find_w_col (graph)
        if (length (w_col) == 0) # sc ensures this never happens, so not covered
            w_col <- d_col # nocov

        fr_col <- find_fr_id_col (graph)
        to_col <- find_to_id_col (graph)

        xfr <- yfr <- xto <- yto <- NA
        # TODO: Modify for other complex but non-spatial types of graph
        if (is_graph_spatial (graph))
        {
            spcols <- find_spatial_cols (graph)
            graph <- tbl_to_df (graph)

            fr_is_num <- vapply (spcols$fr_col, function (i)
                                 is.numeric (graph [[i]]), logical (1))
            to_is_num <- vapply (spcols$to_col, function (i)
                                 is.numeric (graph [[i]]), logical (1))
            if (!(all (fr_is_num) & all (to_is_num)))
                stop (paste0 ("graph appears to have non-numeric ",
                              "longitudes and latitudes"))

            xfr <- spcols$fr_col [1]
            yfr <- spcols$fr_col [2]
            xto <- spcols$to_col [1]
            yto <- spcols$to_col [2]
        } else
        {
            if (length (fr_col) != 1 & length (to_col) != 1)
                stop ("Unable to determine from and to columns in graph") # nocov
        }
    }

    time_col <- grep ("time", nms)
    if (length (time_col) != 1)
    {
        time_col <- grep ("time$", nms)
        if (length (time_col) != 1)
            time_col <- NA
    }
    timew_col <- grep ("time_w|timew|tw", nms)
    if (length (timew_col) != 1)
    {
        timew_col <- grep ("time_w|timew|^tw", nms)
        if (length (timew_col) != 1)
            timew_col <- NA
    }

    ret <- c (edge_id, fr_col, to_col, d_col, w_col, time_col, timew_col,
              xfr, yfr, xto, yto, component)
    names (ret) <- c ("edge_id", "from", "to", "d", "d_weighted", "time",
                      "time_weighted", "xfr", "yfr", "xto", "yto", "component")
    class (ret) <- c (class (ret), "graph_columns")

    # This is passed to many C++ routines, in which case it needs to be
    # converted to a vector (`do.call (c, gr_cols)`), and the R-style 1-indexeso
    # need to be converted to equivalent 0-indexed forms
    return (as.list (ret))
}

#' convert_graph
#'
#' Convert graph to a standard form suitable for submission to C++ routines
#' @noRd
convert_graph <- function (graph, gr_cols)
{
    keep_cols <- c ("edge_id", "from", "to", "d", "d_weighted",
                    "time", "time_weighted")
    index <- do.call (c, gr_cols [keep_cols])
    index <- index [!is.na (index)]
    graph <- graph [, index]
    names (graph) <- names (index)

    if ("edge_id" %in% names (graph))
        graph$edge_id <- convert_to_char (graph$edge_id)
    graph$from <- convert_to_char (graph$from)
    graph$to <- convert_to_char (graph$to)

    if (!"time_weighted" %in% names (graph))
        graph$time_weighted <- graph$time

    return (graph)
}

convert_to_char <- function (x)
{
    if (!is.character (x)) x <- paste0 (x)
    return (x)
}

tbl_to_df <- function (graph)
{
    if (methods::is (graph, "tbl"))
    {
        classes <- class (graph) [!grepl ("tbl", class (graph))]
        graph <- as.data.frame (graph)
        class (graph) <- classes
    }
    return (graph)
}

#' make_vert_map
#'
#' Map unique vertex names to sequential numbers in matrix
#' @noRd
make_vert_map <- function (graph, gr_cols, xy = FALSE)
{
    # gr_cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    verts <- c (paste0 (graph [[gr_cols$from]]), paste0 (graph [[gr_cols$to]]))
    indx <- which (!duplicated (verts))
    if (!xy)
    {
        # Note id has to be 0-indexed:
        res <- data.frame (vert = paste0 (verts [indx]),
                           id = seq (indx) - 1,
                           stringsAsFactors = FALSE)
    } else
    {
        verts_x <- c (graph [[gr_cols$xfr]], graph [[gr_cols$xto]])
        verts_y <- c (graph [[gr_cols$yfr]], graph [[gr_cols$yto]])
        res <- data.frame (vert = paste0 (verts [indx]),
                           id = seq (indx) - 1,
                           x = verts_x [indx],
                           y = verts_y [indx],
                           stringsAsFactors = FALSE)
    }
    return (res)
}

#' get_index_id_cols
#'
#' Get an index of `pts` matching `vert_map`, as well as the
#' corresonding names of those `pts`
#'
#' @return list of `index`, which is 0-based for C++, and corresponding
#' `id` values.
#' @noRd
get_index_id_cols <- function (graph, gr_cols, vert_map, pts)
{
    index <- -1
    id <- NULL
    if (!missing (pts))
    {
        if (is.integer (pts) & is.vector (pts))
        {
            index <- pts
        } else if (is.character (pts) | is.numeric (pts) |
            is.matrix (pts) | is.data.frame (pts))
        {
            index <- get_pts_index (graph, gr_cols, vert_map, pts)
        } else
            stop ("routing points are of unknown form; must be either ",
                  "character, matrix, or integer")

        if (length (pts == 2) & is.numeric (pts) &
            ((any (grepl ("x", names (pts), ignore.case = TRUE)) &
             any (grepl ("y", names (pts), ignore.case = TRUE))) |
             (any (grepl ("lon", names (pts), ignore.case = TRUE) &
                   (any (grepl ("lat", names (pts), ignore.case = TRUE)))))))
            names (pts) <- NULL

        id <- get_id_cols (pts)
        if (is.null (id))
            id <- vert_map$vert [index] # index is 1-based
    }
    list (index = index, id = id)
}

#' get_pts_index
#'
#' Convert `from` or `to` args of `dodgr_dists` to indices into
#' `vert_map`
#'
#' @param graph A dodgr graph
#' @param vert_map Two-column `data.frame` of unique vertices and
#' corresponding IDs, obtained from `make_vert_map`
#' @param gr_cols Returned from `dodgr_graph_cols()`
#' @param pts Either a vector of names, or a matrix or `data.frame` of
#' arbitrary geographical coordinates for which to get index into vertices of
#' graph.
#'
#' @noRd
get_pts_index <- function (graph, gr_cols, vert_map, pts)
{
    if (!(is.matrix (pts) | is.data.frame (pts)))
    {
        if (!is.numeric (pts))
            pts <- matrix (pts, ncol = 1)
        else
        {
            nms <- names (pts)
            if (is.null (nms))
                nms <- c ("x", "y")
            pts <- matrix (pts, nrow = 1) # vector of (x,y) vals
            colnames (pts) <- nms
        }
    }

    if (ncol (pts) == 1)
    {
        pts <- pts [, 1]
        if (!is.numeric (pts))
        {
            indx <- match (pts, vert_map$vert)
            if (any (is.na (indx)))
                stop (paste0 ("from/to are not numeric yet can not be",
                              " matched onto graph vertices"))
            pts <- indx
        }
        if (any (pts < 1 | pts > nrow (vert_map)))
            stop (paste0 ("points exceed numbers of vertices"))
    } else
    {
        nms <- names (pts)
        if (is.null (nms))
            nms <- colnames (pts)
        ix <- which (grepl ("x", nms, ignore.case = TRUE) |
                     grepl ("lon", nms, ignore.case = TRUE))
        iy <- which (grepl ("y", nms, ignore.case = TRUE) |
                     grepl ("lat", nms, ignore.case = TRUE))
        if (length (ix) != 1 | length (iy) != 1)
            stop (paste0 ("Unable to determine geographical ",
                          "coordinates in from/to"))

        index <- match (c ("xfr", "yfr", "xto", "yto"), names (gr_cols))
        if (any (is.na (gr_cols [index])))
            stop (paste0 ("Cannot determine geographical coordinates ",
                          "against which to match pts"))

        if (is.data.frame (pts))
        {
            names (pts) [ix] <- "x"
            names (pts) [iy] <- "y"
        } else
        {
            colnames (pts) [ix] <- "x"
            colnames (pts) [iy] <- "y"
        }

        # Result of rcpp_points_index is 0-indexed for C++
        pts <- rcpp_points_index_par (dodgr::dodgr_vertices (graph), pts) + 1
        # xy has same order as vert_map
    }

    pts
}

#' get_id_cols
#'
#' Get the ID columns or rownames from a matrix or data.frame of from or to
#' points
#'
#' @param pts The `from` or `to` args passed to `dodgr_dists`
#' @return Character vector of names of points, if they exist in `pts`
#' @noRd
get_id_cols <- function (pts)
{
    ids <- NULL
    if (any (grepl ("id", colnames (pts), ignore.case = TRUE)))
    {
        nmc <- which (grepl ("id", colnames (pts)))
        if (methods::is (pts, "data.frame"))
            ids <- pts [[nmc]]
        else if (is.matrix (pts))
            ids <- pts [, nmc, drop = TRUE]
    } else if (is.vector (pts) & !is.null (names (pts)))
        ids <- names (pts)
    else if (!is.null (rownames (pts)))
        ids <- rownames (pts)
    return (ids)
}

