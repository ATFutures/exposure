context("exposure")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

source ("../sc-conversion-fns.R")

test_that("exposure", {
              h <- sf_to_sc (dodgr::hampi)
              net <- exposure_centrality (h)
              expect_is (net, "dodgr_streetnet_sc")
              expect_true ("centrality" %in% names (net))
              expect_true ("exposure" %in% names (net))
})
