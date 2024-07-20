test_that("Utils check - similarity and co", {
  signatures <- list(
    "signature1" = c("gene1", "gene2", "gene3"),
    "signature2" = c("gene2", "gene3", "gene4"),
    "signature3" = c("gene1", "gene5")
  )
  similarity <- signatures_similarity(signatures, metric_name = "jaccard")

  expect_true(is(similarity, "data.frame"))
  expect_true(all(dim(similarity) == c(3, 3)))
  expect_true(all(similarity <= 1))

  sim_perc <- signatures_similarity(signatures, metric_name = "percentage")
  expect_true(all(sim_perc <= 100))
})
