test_that("shared input normalization preserves numeric matrix inputs", {
  x <- matrix(seq_len(20) / 10, nrow = 5L, ncol = 4L)
  colnames(x) <- paste0("x", seq_len(ncol(x)))

  expect_identical(.selsurv_as_design_matrix(x), x)
  expect_identical(as_numeric_matrix(x), x)
})


test_that("categorical coding is treatment coded and invariant to column order", {
  d1 <- data.frame(
    first = factor(c("A", "B", "C", "A", "B", "C")),
    second = factor(c("u", "v", "u", "v", "u", "v")),
    flag = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
    z = seq_len(6)
  )

  d2 <- d1[c("second", "first", "flag", "z")]

  m1 <- .selsurv_as_design_matrix(d1)
  m2 <- .selsurv_as_design_matrix(d2)

  expect_identical(
    colnames(m1),
    c("firstB", "firstC", "secondv", "flagTRUE", "z")
  )

  expect_setequal(colnames(m1), colnames(m2))

  expect_equal(
    m1,
    m2[, colnames(m1), drop = FALSE],
    tolerance = 0,
    check.attributes = FALSE
  )

  expected <- cbind(
    firstB = as.numeric(d1$first == "B"),
    firstC = as.numeric(d1$first == "C"),
    secondv = as.numeric(d1$second == "v"),
    flagTRUE = as.numeric(d1$flag),
    z = as.numeric(d1$z)
  )

  expect_equal(
    unname(m1),
    unname(expected),
    tolerance = 0
  )
})


test_that("ordered and character predictors use ordinary treatment coding", {
  dat <- data.frame(
    binary = factor(
      c("control", "treated", "control", "treated"),
      levels = c("control", "treated")
    ),
    logical = c(FALSE, TRUE, FALSE, TRUE),
    ordered = ordered(
      c("low", "medium", "high", "medium"),
      levels = c("low", "medium", "high")
    ),
    character = c("x", "y", "z", "x"),
    numeric = c(1.2, 2.3, 3.4, 4.5),
    stringsAsFactors = FALSE
  )

  got <- .selsurv_as_design_matrix(dat)

  expect_identical(
    colnames(got),
    c(
      "binarytreated",
      "logicalTRUE",
      "orderedmedium",
      "orderedhigh",
      "charactery",
      "characterz",
      "numeric"
    )
  )

  expect_false(any(grepl("\\.(L|Q|C)$", colnames(got))))

  expect_equal(
    unname(got[, "orderedmedium"]),
    as.numeric(dat$ordered == "medium"),
    tolerance = 0
  )

  expect_equal(
    unname(got[, "orderedhigh"]),
    as.numeric(dat$ordered == "high"),
    tolerance = 0
  )
})


test_that("explicit factor level order determines the reference level", {
  dat <- data.frame(
    group = factor(
      c("A", "B", "C", "A", "B", "C"),
      levels = c("C", "A", "B")
    )
  )

  got <- .selsurv_as_design_matrix(dat)

  expect_identical(colnames(got), c("groupA", "groupB"))
  expect_equal(
    unname(got[, "groupA"]),
    as.numeric(dat$group == "A")
  )
  expect_equal(
    unname(got[, "groupB"]),
    as.numeric(dat$group == "B")
  )
})


test_that("factor coding does not depend on global contrast options", {
  old <- options("contrasts")
  on.exit(options(old), add = TRUE)

  options(contrasts = c("contr.sum", "contr.poly"))

  dat <- data.frame(
    group = factor(c("A", "B", "C", "A", "B", "C")),
    ordered = ordered(
      c("low", "medium", "high", "low", "medium", "high"),
      levels = c("low", "medium", "high")
    )
  )

  got <- .selsurv_as_design_matrix(dat)

  expect_identical(
    colnames(got),
    c("groupB", "groupC", "orderedmedium", "orderedhigh")
  )
})


test_that("single-level categorical predictors fail explicitly", {
  dat_factor <- data.frame(
    group = factor(rep("A", 6L)),
    z = seq_len(6L)
  )

  expect_error(
    .selsurv_as_design_matrix(dat_factor),
    "at least two observed categorical levels"
  )

  dat_character <- data.frame(
    group = rep("A", 6L),
    z = seq_len(6L),
    stringsAsFactors = FALSE
  )

  expect_error(
    .selsurv_as_design_matrix(dat_character),
    "at least two observed categorical levels"
  )
})


test_that("non-finite covariates are rejected centrally", {
  x_inf <- matrix(c(1, 2, 3, Inf), nrow = 2L)

  expect_error(
    .selsurv_as_design_matrix(x_inf),
    "finite"
  )

  x_na <- data.frame(
    a = c(1, NA_real_),
    group = factor(c("A", "B"))
  )

  expect_error(
    .selsurv_as_design_matrix(x_na),
    "finite"
  )
})


test_that("sparse numeric matrices match dense numeric matrices", {
  skip_if_not_installed("Matrix")

  x <- matrix(seq_len(24) / 10, nrow = 8L, ncol = 3L)
  colnames(x) <- c("x1", "x2", "x3")

  xs <- Matrix::Matrix(x, sparse = TRUE)

  expect_equal(
    .selsurv_as_design_matrix(xs),
    x,
    tolerance = 0
  )

  expect_equal(
    as_numeric_matrix(xs),
    x,
    tolerance = 0
  )
})


test_that("procedure-level bootstrap matches explicit treatment-coded input", {
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")

  set.seed(1402)
  n <- 72L

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), length.out = n)),
    flag = rep(c(FALSE, TRUE), length.out = n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )

  x <- cbind(
    groupB = as.numeric(dat$group == "B"),
    groupC = as.numeric(dat$group == "C"),
    flagTRUE = as.numeric(dat$flag),
    z1 = dat$z1,
    z2 = dat$z2
  )
  storage.mode(x) <- "double"

  lp <- 0.50 * dat$z1 -
    0.35 * dat$z2 +
    0.25 * as.numeric(dat$group == "B")

  t_event <- rexp(n, rate = 0.08 * exp(lp))
  t_cens <- rexp(n, rate = 0.05)
  time <- pmin(t_event, t_cens)
  status <- as.integer(t_event <= t_cens)

  support <- match("z1", colnames(x))
  expect_false(is.na(support))

  common <- list(
    time = time,
    status = status,
    B = 1L,
    methods = "fixed_support",
    seed = 91L,
    nfolds = 5L,
    original_support = support,
    original_lambda = 0.05,
    standardize = FALSE,
    ties = "breslow",
    rng_mode = "independent_streams",
    progress = FALSE
  )

  fit_df <- do.call(
    coxlasso_boot_reselect,
    c(list(x = dat), common)
  )

  fit_mm <- do.call(
    coxlasso_boot_reselect,
    c(list(x = x), common)
  )

  expect_identical(fit_df$x_colnames, fit_mm$x_colnames)
  expect_equal(
    unname(fit_df$generator$x),
    unname(fit_mm$generator$x),
    tolerance = 0
  )

  expect_equal(
    fit_df$original$refit$beta,
    fit_mm$original$refit$beta,
    tolerance = 1e-12
  )

  expect_equal(
    fit_df$beta$fixed_support,
    fit_mm$beta$fixed_support,
    tolerance = 1e-12
  )

  expect_equal(
    fit_df$se$fixed_support,
    fit_mm$se$fixed_support,
    tolerance = 1e-12
  )
})


test_that("numeric 0/1 covariates remain numeric rather than categorical", {
  dat <- data.frame(
    binary = c(0, 1, 0, 1),
    z = c(1.2, 2.3, 3.4, 4.5)
  )

  got <- .selsurv_as_design_matrix(dat)

  expect_identical(colnames(got), c("binary", "z"))
  expect_equal(unname(got[, "binary"]), c(0, 1, 0, 1))
})


test_that("unused factor levels are dropped before treatment coding", {
  dat <- data.frame(
    group = factor(
      c("A", "B", "A", "B"),
      levels = c("A", "B", "C")
    ),
    z = 1:4
  )

  got <- .selsurv_as_design_matrix(dat)

  expect_identical(colnames(got), c("groupB", "z"))
  expect_false("groupC" %in% colnames(got))
})


test_that("constant logical predictors fail explicitly", {
  dat_false <- data.frame(
    flag = rep(FALSE, 5L),
    z = seq_len(5L)
  )

  dat_true <- data.frame(
    flag = rep(TRUE, 5L),
    z = seq_len(5L)
  )

  expect_error(
    .selsurv_as_design_matrix(dat_false),
    "at least two observed categorical levels"
  )

  expect_error(
    .selsurv_as_design_matrix(dat_true),
    "at least two observed categorical levels"
  )
})


test_that("unsupported data-frame column classes fail explicitly", {
  dat_date <- data.frame(
    day = as.Date("2026-01-01") + 0:3,
    z = 1:4
  )

  expect_error(
    .selsurv_as_design_matrix(dat_date),
    "unsupported class"
  )

  dat_list <- data.frame(z = 1:4)
  dat_list$bad <- I(as.list(1:4))

  expect_error(
    .selsurv_as_design_matrix(dat_list),
    "unsupported class"
  )
})


test_that("invalid input and generated design names fail explicitly", {
  dup <- data.frame(a = 1:4, b = 5:8)
  names(dup) <- c("x", "x")

  expect_error(
    .selsurv_as_design_matrix(dup),
    "must be unique"
  )

  empty <- data.frame(a = 1:4, b = 5:8)
  names(empty) <- c("", "b")

  expect_error(
    .selsurv_as_design_matrix(empty),
    "non-missing and non-empty"
  )

  tick <- data.frame(a = 1:4)
  names(tick) <- "bad`name"

  expect_error(
    .selsurv_as_design_matrix(tick),
    "backticks"
  )

  reserved <- data.frame(.time = 1:4, z = 5:8)

  expect_error(
    .selsurv_as_design_matrix(reserved),
    "reserved"
  )

  collision <- data.frame(
    group = factor(c("A", "B", "A", "B")),
    groupB = 1:4,
    check.names = FALSE
  )

  expect_error(
    .selsurv_as_design_matrix(collision),
    "must be unique"
  )
})


test_that("empty covariate inputs fail explicitly", {
  expect_error(
    .selsurv_as_design_matrix(matrix(numeric(0), nrow = 0L, ncol = 1L)),
    "at least one row"
  )

  expect_error(
    .selsurv_as_design_matrix(data.frame(row.names = 1:3)),
    "at least one covariate"
  )
})


test_that("coxlasso_refit rejects a changed design-coordinate system", {
  x <- matrix(
    seq_len(20) / 10,
    nrow = 10L,
    ncol = 2L,
    dimnames = list(NULL, c("x1", "x2"))
  )

  fit <- list(
    X = x,
    support = 1L,
    time = seq_len(10L),
    status = rep(c(1L, 0L), 5L)
  )

  expect_error(
    coxlasso_refit(
      fit,
      X = x[, c("x2", "x1"), drop = FALSE],
      time = fit$time,
      status = fit$status
    ),
    "same design-matrix columns"
  )
})


test_that("supplied re-selection original must match design coordinates", {
  x <- matrix(
    seq_len(18) / 10,
    nrow = 6L,
    ncol = 3L,
    dimnames = list(NULL, c("x1", "x2", "x3"))
  )

  original <- list(
    selection = list(
      support = 1L,
      lambda_min = 0.1
    ),
    refit = list(
      success = TRUE,
      beta = setNames(c(0.5, 0, 0), c("x2", "x1", "x3"))
    )
  )

  expect_error(
    .validate_original_workflow(original, x),
    "coordinates must match x exactly"
  )
})
