% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mmm.R
\name{make_fixed_formula}
\alias{make_fixed_formula}
\title{Generate fixed formula for stacked data}
\usage{
make_fixed_formula(covars = NULL)
}
\arguments{
\item{covars}{a list of character vectors with the covariates names}
}
\value{
a list with the fixed formulas for the pairwise model fitting stage.
}
\description{
Takes a list of covariates. Each element of the list are the covariates for a pairwise model
}
\examples{
\dontrun{
fixed <- make_fixed_formula(covars = c("time", "sex", "time_failure"))
}
}
