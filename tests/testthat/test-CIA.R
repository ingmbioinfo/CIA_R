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

  vec_labels <- CIA_classify(
    data = mat_logcounts,
    signatures_input = gmt,
    n_cpus = 2
  )

  expect_vector(vec_labels)
  expect_true(length(vec_labels) == 379)
  expect_true(length(levels(vec_labels)) == 20)


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


test_that("Extra content, majority voting", {
  ## with sce...
  sce_aio <- CIA_classify(
    data = sce,
    signatures_input = gmt,
    n_cpus = 2,
    similarity_threshold = 0.1,
    column_name = "CIA_prediction_t=0.1"
  )

  table(sce_aio$`CIA_prediction_t=0.1`)

  sce_aio$quick_cluster <-
    scran::quickCluster(sce_aio, min.size=5, assay.type = "tophat_counts")


  ## with seurat...
  so_aio <- CIA_classify(
    data = so,
    signatures_input = gmt,
    seurat_assay = "endogenous",
    seurat_layer = "data",
    n_cpus = 2,
    similarity_threshold = 0.1,
    column_name = "CIA_prediction_t=0.1"
  )
  # so_aio <- Seurat::FindVariableFeatures(so_aio)
  # so_aio <- Seurat::ScaleData(so_aio)
  # so_aio <- Seurat::RunPCA(so_aio)
  so_aio <- FindNeighbors(so_aio, reduction = "PCA", dims = 4)
  so_aio_majority <- celltypist_majority_vote(
    so_aio,
    graph_name = "endogenous_snn",
    classification_obs = "CIA_prediction_t=0.1",
    min_prop = 0.5
  )
  table(so_aio$`CIA_prediction_t=0.1`)
  table(so_aio_majority$`CIA_prediction_t=0.1_majority_voting`)

  expect_true(is.factor(so_aio_majority$`CIA_prediction_t=0.1_majority_voting`))

  # sce_aio_majority <- celltypist_majority_vote(
  #   sce_aio,
  #   classification_obs = "CIA_prediction_t=0.1",
  #   groups_obs = "quick_cluster",
  #   min_prop = 0.5
  # )
  ## TODO: this is to be included/tested properly

})


test_that("Checking performance metrics and associated functionality", {
  sce_aio <- CIA_classify(
    data = sce,
    signatures_input = gmt,
    n_cpus = 2,
    similarity_threshold = 0.1,
    column_name = "CIA_prediction_t=0.1"
  )

  sce_aio <- CIA_classify(
    data = sce_aio,
    signatures_input = gmt,
    n_cpus = 2,
    column_name = "CIA_prediction"
  )

  so_aio <- CIA_classify(
    data = so,
    signatures_input = gmt,
    seurat_assay = "endogenous",
    seurat_layer = "data",
    n_cpus = 2,
    similarity_threshold = 0.1,
    column_name = "CIA_prediction_t=0.1"
  )
  so_aio <- CIA_classify(
    data = so_aio,
    signatures_input = gmt,
    seurat_assay = "endogenous",
    seurat_layer = "data",
    n_cpus = 2,
    column_name = "CIA_prediction"
  )
  so_aio_withallinfo <- score_all_signatures(
    data = so_aio,
    signatures_input = gmt,
    seurat_assay = "endogenous",
    seurat_layer = "data",
    n_cpus = 2
  )

  ccm1 <- compute_classification_metrics(
    cells_info = colData(sce_aio),
    classification_cols = c("CIA_prediction"),
    ref_labels = "CIA_prediction"
  )

  ccm2 <- compute_classification_metrics(
    cells_info = colData(sce_aio),
    classification_cols = c("CIA_prediction", "CIA_prediction_t=0.1"),
    ref_labels = "CIA_prediction",
    unassigned_label = "Unassigned"
  )

  expect_type(ccm1, "double")
  expect_type(ccm2, "double")



  ccm3 <- compute_classification_metrics(
    cells_info = so_aio@meta.data,
    classification_cols = c("CIA_prediction", "CIA_prediction_t=0.1"),
    ref_labels = "CIA_prediction",
    unassigned_label = "Unassigned"
  )
  expect_type(ccm3, "double")


  gccm1 <- grouped_classification_metrics(
    cells_info = colData(sce_aio),
    classification_col= "CIA_prediction",
    ref_labels = "CIA_prediction",
    unassigned_label = "Unassigned"
  )
  expect_type(gccm1, "list")

  gccm3 <- grouped_classification_metrics(
    cells_info = so_aio@meta.data,
    classification_col= "CIA_prediction",
    ref_labels = "CIA_prediction",
    unassigned_label = "Unassigned"
  )
  expect_type(gccm3, "list")

  # Working on the plots -------------------------------------------------------
  p_gc <- plot_group_composition(
    so_aio@meta.data,
    comp_col = "CIA_prediction", ref_col = "CIA_prediction",
    palette = colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(20)
  )
  expect_true(is(p_gc, "gg"))

  p_gch <- group_composition(so_aio@meta.data,
                             classification_obs = "CIA_prediction",
                             ref_obs = "CIA_prediction")
  expect_true(is(p_gch, "gg"))

  p_gd <- grouped_distributions(
    so_aio_withallinfo@meta.data,
    ref_obs = "CIA_prediction",
    columns_obs = names(gmt))
  expect_true(is(p_gd, "gg"))
})

# TODO: rename the param column_name? cia_name?

# TODO: option to ALSO store the CIA scores?

