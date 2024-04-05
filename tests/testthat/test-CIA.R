test_that("Loading signature files", {
  expect_error(
    load_signatures("file_not_found")
  )

  expect_error(
    load_signatures(42)
  )
})


test_that("Individual scores computations", {
  # individual - SCE --------------------------------------------------------
  scores <- score_signature(
    data = sce,
    sce_assay = "logcounts",
    geneset = gmt$`L5 ET`
  )
  colData(sce)$score_L5et <- scores
  expect_vector(scores)

  nullscores <- score_signature(
    data = sce,
    sce_assay = "logcounts",
    geneset = ""
  )
  expect_true(all(nullscores == 0))

  expect_error(
    score_signature(
      data = "my_data",
      sce_assay = "logcounts",
      geneset = gmt$`L5 ET`
    )
  )

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

  # individual - Seurat -----------------------------------------------------
  ssc <- score_signature(
    data = so,
    seurat_assay = "endogenous",
    seurat_layer = "data",
    geneset = gmt$`L5 ET`
  )
  expect_vector(ssc)
  expect_true(identical(scores, ssc))

  # individual - matrix -----------------------------------------------------
  dim(mat_logcounts)
  sc2 <- score_signature(
    data = mat_logcounts,
    geneset = gmt$`L5 ET`
  )

  expect_vector(sc2)
  expect_true(identical(scores, sc2))

  # Empty matrix
  expect_error(
    score_signature(
      data = matrix(),
      geneset = gmt$`L5 ET`
    )
  )

  expect_error(
    score_signature(
      data = as.matrix(numeric(0)),
      geneset = gmt$`L5 ET`
    )
  )

  expect_error(
    score_signature(
      data = Matrix(),
      geneset = gmt$`L5 ET`
    )
  )

  # Matrix of chars...
  expect_error(
    score_signature(
      data = matrix(data = letters),
      geneset = gmt$`L5 ET`
    )
  )

})



test_that("All scores at once computations", {

  # all signatures - SCE ----------------------------------------------------
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


  # all signatures - Seurat -------------------------------------------------
  so_cia <- score_all_signatures(
    data = so,
    signatures_input = gmt,
    seurat_assay = "endogenous",
    seurat_layer = "data",
    return_score = FALSE,
    score_mode = "scaled",
    n_cpus = 2
  )



  # all signatures - matrix -------------------------------------------------
  # TODO

  allsigs_cia <- score_all_signatures(
    data = mat_logcounts,
    signatures_input = gmt,
    return_score = FALSE,
    score_mode = "scaled",
    n_cpus = 2
  )

  expect_true(is.matrix(allsigs_cia))
  expect_identical(
    dim(allsigs_cia),
    c(ncol(mat_logcounts), length(gmt))
  )

  allsigs_cia_logged <- score_all_signatures(
    data = mat_logcounts,
    signatures_input = gmt,
    return_score = FALSE,
    score_mode = "log10",
    n_cpus = 2
  )
  expect_true(is.matrix(allsigs_cia_logged))
  expect_identical(
    dim(allsigs_cia_logged),
    c(ncol(mat_logcounts), length(gmt))
  )

  # Running the all in one function...

  sce_aio <- CIA_classify(
    data = sce_cia,
    signatures_input = gmt,
    n_cpus = 2,
    similarity_threshold = 0.1,
    column_name = "CIA_prediction_t=0.1"
  )


  so_aio <- CIA_classify(
    data = so_cia,
    signatures_input = gmt,
    seurat_assay = "endogenous",
    seurat_layer = "data",
    n_cpus = 2,
    similarity_threshold = 0.1,
    column_name = "CIA_prediction_t=0.1"
  )


  expect_error(
    CIA_classify(
      data = 42,
      signatures_input = gmt,
      n_cpus = 2,
      similarity_threshold = 0.1,
      column_name = "CIA_prediction_t=0.1"
    )
  )

  expect_error(
    CIA_classify(
      data = sce_cia,
      signatures_input = "file_not_found",
      n_cpus = 2,
      similarity_threshold = 0.1,
      column_name = "CIA_prediction_t=0.1"
    )
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

