test_that("repermulize_wrapper runs on a small simulated example", {
  set.seed(1)

  tree <- ape::rtree(8)

  pheno <- stats::rnorm(8)
  names(pheno) <- tree$tip.label

  genes <- matrix(
    sample(c(0, 1), size = 5 * 8, replace = TRUE),
    nrow = 5,
    dimnames = list(paste0("gene", 1:5), tree$tip.label)
  )

  res <- repermulize_wrapper(
    real_pheno = pheno,
    real_genes = genes,
    real_tree = tree,
    n = 20,
    use_futures = FALSE,
    verbose = FALSE,
    regression_method = "lm"
  )

  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), nrow(genes))
  expect_true(all(c("Estimate", "p.value") %in% rownames(res)))
})
