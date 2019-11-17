context("exposure")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

source ("../sc-conversion-fns.R")

test_that("exposure", {
              net <- sf_to_sc (dodgr::hampi)
              net_v <- dodgr::weight_streetnet (net, wt_profile = "motorcar")
              net_v <- exposure_centrality (net_v)
              expect_is (net_v, "dodgr_streetnet_sc")
              expect_true ("centrality" %in% names (net_v))
              net_p <- dodgr::weight_streetnet (net, wt_profile = "foot")
              net_p <- exposure (net_p, net_v)
              expect_is (net_p, "dodgr_streetnet_sc")
              expect_true ("exposure" %in% names (net_p))
})
