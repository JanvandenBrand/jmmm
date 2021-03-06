% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_mmm.R
\name{get_likelihood_profiles_nofail}
\alias{get_likelihood_profiles_nofail}
\title{Get likelihood profiles for non failures}
\usage{
get_likelihood_profiles_nofail(likelihood_nofail, landmark, horizon, interval)
}
\arguments{
\item{likelihood_nofail}{a list of likelihoods for each subject at the observed time points returned by get_likelihood_nofail().}

\item{landmark}{the prediction landmark.}

\item{horizon}{the prediction horizon.}

\item{interval}{the intervals relative to the horizon (e.g. 1/12 if monthly intervals if horizon is in years)}
}
\value{
A list of length nrow(unique(id)) with likelihood profiles
}
\description{
Generate a step function of the likelihood at interval times.
}
