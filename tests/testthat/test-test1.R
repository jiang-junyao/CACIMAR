test_that("test network", {
  O1 <- OrthG_Mm_Zf
  load(system.file("extdata", "network_example.rda", package = "CACIMAR"))
  n1=Identify_ConservedNetworks(O1,mmNetwork,zfNetwork,'mm','zf')
  expect_equal(length(n1),2)
})

test_that("test CACIMAR", {
  data("pbmc_small")
  all.markers <- Identify_Markers(pbmc_small)
  all.markers2 <- Format_Markers_Frac(all.markers)
  load(system.file("extdata", "zf_mm_markers.rda", package = "CACIMAR"))
  ConservedMarker <- Identify_ConservedMarkers(OrthG_Mm_Zf,Mm_marker,Zf_marker,
                                               Species_name1 = 'mm',Species_name2 = 'zf')
  expression <- Identify_ConservedCellTypes(OrthG_Mm_Zf,Zf_marker[1:100,],Mm_marker[1:100,],'zf','mm')
})
