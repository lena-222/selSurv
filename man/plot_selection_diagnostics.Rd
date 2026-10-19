\name{plot_selection_diagnostics}
\alias{plot_selection_diagnostics}
\title{Plot selection diagnostics across bootstrap workflows}
\usage{
plot_selection_diagnostics(
  diagnostics,
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
  legend_position = c(
    "top_panel", "inside", "none", "side"
  ),
  colors = NULL,
  pch = c(16, 17, 15, 18),
  cex_points = 1.1,
  cex_names = 0.85,
  cex_legend = 0.85,
  las = 2,
  ...
)
}
\arguments{
\item{diagnostics}{A data frame returned by \code{compare_selection_diagnostics()}.}
\item{metric}{Which quantity to plot: \code{"reselect"} or \code{"zero"}.}
\item{plot_type}{Plot type. \code{"dot"} gives a horizontal dot plot; \code{"bar"} gives grouped bars.}
\item{workflows}{Character vector selecting workflows to show.}
\item{sort_by}{Optional column name used to order terms.}
\item{decreasing}{Logical; should sorting be decreasing?}
\item{palette}{Color palette. \code{"color"} gives a colorblind-friendly palette; \code{"grey"} gives neutral greys.}
\item{ylim}{Optional y-axis limits for bar plots.}
\item{xlim}{Optional x-axis limits for dot plots.}
\item{main}{Optional plot title.}
\item{xlab}{Optional x-axis label for dot plots.}
\item{ylab}{Optional y-axis label for bar plots.}
\item{legend}{Logical; add a legend?}
\item{legend_position}{Legend placement. \code{"top_panel"} uses a separate legend panel above the plot and is robust in the RStudio zoom window.}
\item{colors}{Optional vector of colors. Overrides \code{palette}.}
\item{pch}{Point symbols for dot plots.}
\item{cex_points}{Point size for dot plots.}
\item{cex_names}{Axis-label size for term labels.}
\item{cex_legend}{Legend text size.}
\item{las}{Axis label orientation for bar plots.}
\item{...}{Additional graphical arguments passed to plotting functions.}
}
\value{Invisibly returns the matrix used for plotting.}
\description{
Displays re-selection or zero frequencies as dot plots or grouped bars.

Supply a combined selection-diagnostics table (see
\code{compare_selection_diagnostics()}).
}
