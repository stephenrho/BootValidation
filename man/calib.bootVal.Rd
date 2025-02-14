% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallel_vboot.R
\name{calib.bootVal}
\alias{calib.bootVal}
\title{Plot non-linear (lowess) calibration curve}
\usage{
calib.bootVal(x)
}
\arguments{
\item{x}{A bootVal object (must have been run with save_ps = TRUE)}
}
\value{
A plot of apparent and optimism corrected calibration curves

a list of original lowess smooth and optimisms
}
\description{
Plots non-linear calibration curve using lowess smoother
   Also plots optimism corrected calibration curve using method described in Harrell (2015) p. 270
   modified from `rms::calibrate`.
}
\references{
Harrell Jr, F. E. (2015). Regression modeling strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. Springer.
}
