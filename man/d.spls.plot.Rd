% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d.spls.plot.R
\name{d.spls.plot}
\alias{d.spls.plot}
\title{Plots the coefficient curve of a Dual-SPLS regression}
\usage{
d.spls.plot(mod.dspls,ncomp)
}
\arguments{
\item{mod.dspls}{is a fitted Dual-SPLS object.}

\item{ncomp}{a positive integer or a numeric vector of the number of Dual-SPLS components to consider.}
}
\value{
no return value
}
\description{
The function \code{dual.spls.plot} provides the regression coefficient curves of a Dual-SPLS model for a specified number of components
and the mean of the original data plot.
}
\details{
The plots allow the visualization of the results by comparing the mean of the original data to the coefficients regression in the
case of a Dual-SPLS regression. The plots provided correspond to the Dual-SPLS coefficients for each \code{ncomp} desired.
}
\examples{
### load dual.spls library
library(dual.spls)
### constructing the simulated example
n <- 100
p <- 50
nondes <- 20
sigmaondes <- 0.5
data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)

X <- data$X
y <- data$y

#fitting the Dual-SPLS lasso model

mod.dspls.lasso <- d.spls.lasso(X=X,y=y,ncp=10,ppnu=0.9,verbose=TRUE)

ncomp=c(5,6,7)
d.spls.plot(mod.dspls.lasso,ncomp)
}
\author{
Louna Alsouki François Wahl
}
