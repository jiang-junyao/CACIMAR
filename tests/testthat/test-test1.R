test_that("test CACIMAR", {
  O1 <- OrthG_Mm_Zf
  load(system.file("extdata", "network_example.rda", package = "CACIMAR"))
  n1=Identify_ConservedNetworks(O1,mmNetwork,zfNetwork,'mm','zf')
  expect_equal(length(n1),2)
})
