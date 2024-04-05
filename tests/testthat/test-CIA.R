test_that("Individual scores computations", {
  scores <- score_signature(
    data = sce,
    sce_assay = "logcounts",
    geneset = gmt$`L5 ET`
  )
  colData(sce)$score_L5et <- scores


  expect_vector(scores)

  expect_error(
    score_signature(
      data = sce,
      sce_assay = "wrongassay",
      geneset = gmt$`L5 ET`
    )
  )

  expect_error(
    score_signature(
      data = sce$NREADS,
      sce_assay = "logcounts",
      geneset = gmt$pippo
    )
  )

  expect_error(
    score_signature(
      data = sce,
      sce_assay = "logcounts",
      geneset = gmt$pippo
    )
  )


})



test_that("All scores at once computations", {
  sce_cia <- score_all_signatures(
    data = sce,
    signatures_input = gmt,
    return_score = FALSE,
    score_mode = "scaled",
    n_cpus = 2
  )

  expect_true(is(sce_cia, "SingleCellExperiment"))

  expect_error(
    score_all_signatures(
      data = colData(sce_cia),
      signatures_input = gmt,
      return_score = FALSE,
      score_mode = "scaled",
      n_cpus = 2
    )
  )

  expect_error(
    score_all_signatures(
      data = sce,
      signatures_input = gmt$Astro,
      return_score = FALSE,
      score_mode = "scaled",
      n_cpus = 2
    )
  )

  expect_error(
    score_all_signatures(
      data = sce,
      signatures_input = gmt,
      return_score = "yes",
      score_mode = "scaled",
      n_cpus = 2
    )
  )

  expect_error(
    score_all_signatures(
      data = sce,
      signatures_input = gmt,
      return_score = FALSE,
      score_mode = 1,
      n_cpus = 2
    )
  )

  expect_error(
    score_all_signatures(
      data = sce,
      signatures_input = gmt,
      return_score = FALSE,
      score_mode = "scaled",
      n_cpus = "parallel"
    )
  )



  # Running the all in one function...

  sce_aio <- CIA_classify(
    data = sce_cia,
    signatures_input = gmt,
    n_cpus = 2,
    similarity_threshold = 0.1,
    column_name = "CIA_prediction_t=0.1"
  )


  plotTSNE(sce_aio, colour_by = "CIA_prediction_t=0.1")

})





# TODO: rename the param column_name? cia_name?

# TODO: option to ALSO store the CIA scores?

# compute_classification_metrics(
#   cells_info = colData(sce_aio),
#   classification_cols = c("CIA_prediction_t=0.1"),
#   ref_labels = "Primary.Type"
# )
#
# # removing unassigned from metrics calculation
# compute_classification_metrics(
#   cells_info = colData(sce_aio),
#   classification_cols = c("CIA_prediction", "CIA_prediction_t=0.1"),
#   ref_labels = "Primary.Type",
#   unassigned_label = "Unassigned"
# )

