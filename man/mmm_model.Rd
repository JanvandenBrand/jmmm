% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mmm.R
\name{mmm_model}
\alias{mmm_model}
\title{Fit a multivariate mixed model using the pairwise fitting approach}
\usage{
mmm_model(
  fixed,
  random,
  id,
  data,
  stacked_data,
  pairs,
  model_families,
  iter_EM = 100,
  iter_qN_outer = 10,
  ncores = NULL,
  ...
)
}
\arguments{
\item{fixed}{character string with the fixed model part. This is passed to as.formula().}

\item{random}{character string with the random model part. This is passed to as.formula().}

\item{id}{The name of the column with subject ids in stacked_data.}

\item{data}{the original data (in wide, unstacked format).}

\item{stacked_data}{a list with the stacked data returned by the stack_data() function.}

\item{pairs}{a character matrix with the pairs returned by the make_pairs() function.}

\item{model_families}{a list with family names and indicators returned by the test_input_datatypes() function.}

\item{iter_EM}{(default = 100) EM iterations passed to GLMMadaptive::mixed_model().}

\item{iter_qN_outer}{(default = 10) outer Quasi Newton iterations passed to GLMMadaptive::mixed_model().}

\item{ncores}{number of workers for multisession parallel plan.}

\item{...}{additional arguments passed to GLMMadaptive::mixed_model().}
}
\value{
A list of 3 elements including:
\itemize{
\item{estimates}{A data.frame with the parameter estimates and standard errors}
\item{vcov}{The Variance-Covariance matrix for the random effects}
\item{corr}{A correlation matrix for the random effects}
}
}
\description{
Takes the stacked data and the user defined formulas for the fixed and random effects. Leverages the future.apply package to parallelize along pairs.
Using a different number of workers than pairs may result in problems during the optimization process. The output of the pairwise model fitting is combined
and any parameters that have been estimated multiple times are averaged. The procedure takes the individual level contributions to the Maximum Likelihood in
order to estimate the standard errors for the parameter estimates with a robust sandwich estimator.
}
\examples{
\dontrun{
# enable parallel processing on a Windows machine
future::plan(multisession, workers=2)

# Fit the model
mmm_model(fixed = fixed,
random = random,
id = "id",
data = df,
stacked_data = df_stacked,
pairs = pairs,
model_families = model_info,
iter_EM = 100,
iter_qN_outer = 30,
nAGQ = 7)

# Reset the parallel plan back to its default.
future::plan("default")
}
}
