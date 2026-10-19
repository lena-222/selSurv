test_that("compare_selection_diagnostics and plot_selection_diagnostics work", {
  d_fixed <- data.frame(
    term = c("x1", "x2"),
    reselect_frequency = c(1, 0),
    zero_frequency = c(0, 1)
  )
  d_t2 <- data.frame(
    term = c("x1", "x2"),
    reselect_frequency = c(0.9, 0.4),
    zero_frequency = c(0.1, 0.6)
  )

  cmp <- compare_selection_diagnostics(
    fixed = list(diagnostics = d_fixed),
    reselect_fixed_lambda = list(diagnostics = d_t2),
    digits = 2
  )

  expect_true(all(c("fixed_reselect", "t2_reselect") %in% names(cmp)))
  expect_equal(nrow(cmp), 2)

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  mat <- plot_selection_diagnostics(cmp, metric = "reselect", legend = FALSE)
  mat_bar <- plot_selection_diagnostics(cmp, metric = "reselect", plot_type = "bar", legend = FALSE)
  grDevices::dev.off()

  expect_true(is.matrix(mat))
  expect_true(is.matrix(mat_bar))
  expect_true(file.exists(tf))
})
