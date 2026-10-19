#' Plot selection diagnostics across bootstrap workflows
#'
#' Displays re-selection or zero frequencies as dot plots or grouped bars.
#'
#' Supply a combined selection-diagnostics table (see
#' [compare_selection_diagnostics()]).
#'
#' @param diagnostics A data frame returned by [compare_selection_diagnostics()].
#' @param metric Which quantity to plot: `"reselect"` or `"zero"`.
#' @param plot_type Plot type. `"dot"` gives a horizontal dot plot; `"bar"`
#'   gives grouped bars.
#' @param workflows Character vector selecting workflows to show.
#' @param sort_by Optional column name used to order terms.
#' @param decreasing Logical; should sorting be decreasing?
#' @param palette Color palette. `"color"` gives a colorblind-friendly palette;
#'   `"grey"` gives neutral greys.
#' @param ylim Optional y-axis limits for bar plots. Ignored for dot plots.
#' @param xlim Optional x-axis limits for dot plots.
#' @param main Optional plot title.
#' @param xlab Optional x-axis label for dot plots.
#' @param ylab Optional y-axis label for bar plots.
#' @param legend Logical; add a legend?
#' @param legend_position Legend placement. `"top_panel"` uses a separate legend
#'   panel above the plot and is robust in the RStudio zoom window.
#'   `"inside"` places the legend inside the plotting region. `"none"` suppresses
#'   the legend. `"side"` uses a right-side panel.
#' @param colors Optional vector of colors. Overrides `palette`.
#' @param pch Point symbols for dot plots.
#' @param cex_points Point size for dot plots.
#' @param cex_names Axis-label size for term labels.
#' @param cex_legend Legend text size.
#' @param las Axis label orientation for bar plots.
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @return Invisibly returns the matrix used for plotting.
#' @export
plot_selection_diagnostics <- function(diagnostics,
                                       metric = c("reselect", "zero"),
                                       plot_type = c("dot", "bar"),
                                       workflows = NULL,
                                       sort_by = NULL,
                                       decreasing = TRUE,
                                       palette = c("color", "grey"),
                                       ylim = c(0, 1),
                                       xlim = c(0, 1),
                                       main = NULL,
                                       xlab = NULL,
                                       ylab = NULL,
                                       legend = TRUE,
                                       legend_position = c("top_panel", "inside", "none", "side"),
                                       colors = NULL,
                                       pch = c(16, 17, 15, 18),
                                       cex_points = 1.1,
                                       cex_names = 0.85,
                                       cex_legend = 0.85,
                                       las = 2,
                                       ...) {
  metric <- match.arg(metric)
  plot_type <- match.arg(plot_type)
  palette <- match.arg(palette)
  legend_position <- match.arg(legend_position)

  if (!is.data.frame(diagnostics)) {
    stop("`diagnostics` must be a data frame.", call. = FALSE)
  }
  if (!("term" %in% names(diagnostics))) {
    stop("`diagnostics` must contain a `term` column.", call. = FALSE)
  }

  suffix <- paste0("_", metric)
  value_cols <- grep(paste0(suffix, "$"), names(diagnostics), value = TRUE)
  if (length(value_cols) == 0) {
    stop("No columns ending in `", suffix, "` found.", call. = FALSE)
  }

  if (!is.null(workflows)) {
    wanted_cols <- paste0(workflows, suffix)
    missing <- setdiff(wanted_cols, value_cols)
    if (length(missing) > 0) {
      stop("Requested workflow columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    value_cols <- wanted_cols
  }

  d <- diagnostics
  if (!is.null(sort_by)) {
    if (!(sort_by %in% names(d))) {
      stop("`sort_by` is not a column in `diagnostics`.", call. = FALSE)
    }
    d <- d[order(d[[sort_by]], decreasing = decreasing), , drop = FALSE]
  }

  mat <- t(as.matrix(d[, value_cols, drop = FALSE]))
  storage.mode(mat) <- "numeric"
  colnames(mat) <- d$term
  rownames(mat) <- sub(suffix, "", value_cols, fixed = TRUE)
  rownames(mat) <- sub("^t2$", "Target-2 fixed lambda", rownames(mat))
  rownames(mat) <- sub("^fixed$", "Fixed support", rownames(mat))
  rownames(mat) <- sub("^retuned$", "Retuned", rownames(mat))

  if (is.null(main)) {
    main <- if (metric == "reselect") "Bootstrap re-selection frequencies" else "Bootstrap zero frequencies"
  }
  if (is.null(xlab)) {
    xlab <- if (metric == "reselect") "Re-selection frequency" else "Zero frequency"
  }
  if (is.null(ylab)) {
    ylab <- if (metric == "reselect") "Re-selection frequency" else "Zero frequency"
  }

  if (is.null(colors)) {
    colors <- if (palette == "color") {
      c("#1F77B4", "#D55E00", "#009E73", "#CC79A7", "#7570B3")
    } else {
      grDevices::gray.colors(max(3, nrow(mat)), start = 0.15, end = 0.75)
    }
  }
  if (length(colors) < nrow(mat)) {
    colors <- rep(colors, length.out = nrow(mat))
  }
  if (length(pch) < nrow(mat)) {
    pch <- rep(pch, length.out = nrow(mat))
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    graphics::layout(matrix(1))
  }, add = TRUE)

  use_top_panel <- isTRUE(legend) && legend_position == "top_panel"
  use_side_legend <- isTRUE(legend) && legend_position == "side"

  if (use_top_panel) {
    graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(0.35, 4))
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    if (plot_type == "dot") {
      graphics::legend(
        "center",
        legend = rownames(mat),
        bty = "n",
        cex = cex_legend,
        col = colors[seq_len(nrow(mat))],
        pch = pch[seq_len(nrow(mat))],
        horiz = TRUE,
        x.intersp = 0.7,
        y.intersp = 0.8
      )
    } else {
      graphics::legend(
        "center",
        legend = rownames(mat),
        bty = "n",
        cex = cex_legend,
        fill = colors[seq_len(nrow(mat))],
        horiz = TRUE,
        x.intersp = 0.7,
        y.intersp = 0.8
      )
    }
  } else if (use_side_legend) {
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(3.8, 1.6))
  }

  draw_side_legend <- function() {
    graphics::par(mar = c(2, 0.2, 3, 0.2))
    graphics::plot.new()
    if (plot_type == "dot") {
      graphics::legend(
        "center",
        legend = rownames(mat),
        bty = "n",
        cex = cex_legend,
        col = colors[seq_len(nrow(mat))],
        pch = pch[seq_len(nrow(mat))],
        x.intersp = 0.7,
        y.intersp = 1.1
      )
    } else {
      graphics::legend(
        "center",
        legend = rownames(mat),
        bty = "n",
        cex = cex_legend,
        fill = colors[seq_len(nrow(mat))],
        x.intersp = 0.7,
        y.intersp = 1.1
      )
    }
  }

  if (plot_type == "dot") {
    n_terms <- ncol(mat)
    n_workflows <- nrow(mat)
    y_base <- rev(seq_len(n_terms))
    offsets <- if (n_workflows == 1) 0 else seq(-0.18, 0.18, length.out = n_workflows)

    left_margin <- max(6, min(12, max(nchar(colnames(mat))) * 0.45))
    graphics::par(mar = c(4.2, left_margin, 3.2, 1.2), xpd = FALSE)

    graphics::plot(
      NA,
      xlim = xlim,
      ylim = c(0.5, n_terms + 0.5),
      axes = FALSE,
      xlab = xlab,
      ylab = "",
      main = main,
      ...
    )

    graphics::abline(v = seq(0, 1, by = 0.25), lty = 3, col = "grey88")
    graphics::axis(1, at = seq(0, 1, by = 0.25))
    graphics::axis(2, at = y_base, labels = colnames(mat), las = 1, tick = FALSE, cex.axis = cex_names)
    graphics::box(bty = "l")

    for (i in seq_len(n_workflows)) {
      y_i <- y_base + offsets[i]
      graphics::segments(0, y_i, mat[i, ], y_i, col = "grey82", lwd = 1)
      graphics::points(mat[i, ], y_i, pch = pch[i], col = colors[i], bg = colors[i], cex = cex_points)
    }

    if (legend && legend_position == "inside") {
      graphics::legend(
        "bottomright",
        legend = rownames(mat),
        bty = "o",
        bg = "white",
        box.col = "grey85",
        cex = cex_legend,
        col = colors[seq_len(nrow(mat))],
        pch = pch[seq_len(nrow(mat))]
      )
    }

    if (use_side_legend) draw_side_legend()
    return(invisible(mat))
  }

  # Grouped bar plot
  bottom_margin <- max(6, min(12, ncol(mat) * 0.45))
  graphics::par(mar = c(bottom_margin, 4.3, 3.2, 1.2), xpd = FALSE)

  bp <- graphics::barplot(
    mat,
    beside = TRUE,
    ylim = ylim,
    las = las,
    main = main,
    ylab = ylab,
    col = colors[seq_len(nrow(mat))],
    border = "grey25",
    cex.names = cex_names,
    ...
  )

  graphics::abline(h = seq(0, 1, by = 0.25), lty = 3, col = "grey85")

  if (legend && legend_position == "inside") {
    graphics::legend(
      "topright",
      legend = rownames(mat),
      bty = "o",
      bg = "white",
      box.col = "grey85",
      cex = cex_legend,
      fill = colors[seq_len(nrow(mat))]
    )
  }

  if (use_side_legend) draw_side_legend()
  invisible(mat)
}
