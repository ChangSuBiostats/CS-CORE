% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSCORE_IRLS_cpp.R
\name{set_D}
\alias{set_D}
\title{Set the design matrix for moment-based regressions}
\usage{
set_D(s, D, adjust_setting, covariate_level)
}
\arguments{
\item{s}{A numeric vector of sequencing depths (for mean regression) or squared sequencing depths (for variance and covariance)}

\item{D}{A numeric matrix of intercept and covariates (\code{n x K})}

\item{adjust_setting}{Logical; whether to adjust for covariates}

\item{covariate_level}{A character string indicating whether covariates are assumed to affect
the underlying gene expression levels (\code{"z"}) or the observed counts (\code{"x"}).}
}
\value{
Design matrix (n by K) for moment-based regressions
}
\description{
Set the design matrix for moment-based regressions
}
\keyword{internal}
