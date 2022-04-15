test_that("test CACIMAR", {
  load(system.file("extdata", "network_example.rda", package = "CACIMAR"))
  n1 <- Identify_ConservedNetworks(OrthG_Mm_Zf,mmNetwork,zfNetwork,'mm','zf')
  n2 <- Identify_ConservedNetworks(OrthG_Mm_Zf,zfNetwork,mmNetwork,'zf','mm')
  expect_equal(length(n1),length(n2))
  expect_equal(nrow(n1[[1]]),nrow(n2[[1]]))
})
