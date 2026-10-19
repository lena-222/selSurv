\name{print-methods}
\alias{print.selSurv_fit}
\alias{print.selSurv_boot}
\alias{print.selSurv_debias}
\alias{print.selSurv_split}
\title{Print methods for selSurv objects}
\usage{
\method{print}{selSurv_fit}(x, ...)

\method{print}{selSurv_boot}(x, ...)

\method{print}{selSurv_debias}(x, ...)

\method{print}{selSurv_split}(x, ...)
}
\arguments{
\item{x}{An object returned by a selSurv workflow function.}

\item{...}{Additional arguments, currently ignored.}
}
\description{
Print methods for objects returned by \code{coxlasso_fit()},
\code{coxlasso_boot_parametric()}, \code{coxlasso_debias()}, and
\code{coxlasso_split()}.
}
\value{
The input object \code{x}, invisibly.
}
