\name{compare_selection_diagnostics}
\alias{compare_selection_diagnostics}
\title{Compare selection diagnostics across bootstrap workflows}
\usage{
compare_selection_diagnostics(
  fixed = NULL,
  reselect_fixed_lambda = NULL,
  reselect_retuned = NULL,
  digits = 2
)
}
\arguments{
\item{fixed}{Optional result from \code{coxlasso_boot_parametric(..., workflow = "fixed_support")}.}
\item{reselect_fixed_lambda}{Optional result from \code{coxlasso_boot_parametric(..., workflow = "reselect_fixed_lambda")}.}
\item{reselect_retuned}{Optional result from \code{coxlasso_boot_parametric(..., workflow = "reselect_retuned")}.}
\item{digits}{Number of digits used for rounding. Use \code{NULL} for no rounding.}
}
\value{
A data frame with one row per term and workflow-specific re-selection and zero-frequency columns.
}
\description{
Combines the \code{diagnostics} tables returned by \code{coxlasso_boot_parametric()} for fixed-support,
fixed-lambda re-selection, and retuned re-selection workflows.
}
