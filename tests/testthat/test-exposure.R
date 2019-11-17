context("exposure")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

test_that("exposure", {
    net <- dodgr::weight_streetnet (dodgr::hampi, wt_profile = "foot")
    netc <- dodgr::dodgr_contract_graph (net)
    expect_silent (netc <- exposure_centrality (netc))
    expect_is (netc, "dodgr_streetnet")
    expect_is (netc, "dodgr_contracted")
    expect_true ("centrality" %in% names (netc))
    expect_true (all (is.finite (netc$centrality)))
})
