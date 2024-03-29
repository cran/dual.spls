% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d.spls.pls.R
\name{d.spls.pls}
\alias{d.spls.pls}
\title{Univariate Partial Least Squares regression}
\usage{
d.spls.pls(X,y,ncp,verbose=TRUE)
}
\arguments{
\item{X}{a numeric matrix of predictors values of dimension \code{(n,p)}. Each row represents one observation and each column one predictor variable.}

\item{y}{a numeric vector or a one column matrix of responses. It represents the response variable for each observation.}

\item{ncp}{a positive integer. \code{ncp} is the number of Dual-SPLS components.}

\item{verbose}{a Boolean value indicating whether or not to display the iterations steps. Default value is \code{TRUE}.}
}
\value{
A \code{list} of the following attributes
\item{Xmean}{the mean vector of the predictors matrix \code{X}.}
\item{scores}{the matrix of dimension \code{(n,ncp)} where \code{n} is the number of observations.The \code{scores} represents
the observations in the new component basis computed by the compression step
of the Dual-SPLS.}
\item{loadings}{the matrix of dimension \code{(p,ncp)} that represents the Dual-SPLS components.}
\item{Bhat}{the matrix of dimension \code{(p,ncp)} that regroups the regression coefficients for each component.}
\item{intercept}{the vector of intercept values for each component.}
\item{fitted.values}{the matrix of dimension \code{(n,ncp)} that represents the predicted values of \code{y}}
\item{residuals}{the matrix of dimension \code{(n,ncp)} that represents the residuals corresponding
to the difference between the responses and the fitted values.}
\item{type}{a character specifying the Dual-SPLS norm used. In this case it is \code{ridge}. }
}
\description{
The function \code{d.spls.pls} performs the PLS1 dimensional reduction methodology using Wold's NIPALS algorithm. It is
a Dual-SPLS regression with the norm \deqn{\Omega(w)=||w||_2.}
}
\details{
The resulting solution for \eqn{w} and hence for the coefficients vector, in the PLS regression for one component is
\eqn{w=X^Ty}. In order to compute the next components, a deflation step is performed only considering the parts of
\eqn{X} that are orthogonal to the previous components.
}
\examples{
### load dual.spls library
library(dual.spls)
### constructing the simulated example
oldpar <- par(no.readonly = TRUE)
n <- 100
p <- 50
nondes <- 20
sigmaondes <- 0.5
data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)

X <- data$X
y <- data$y

#fitting the model
mod.dspls <- d.spls.pls(X=X,y=y,ncp=10,verbose=TRUE)

str(mod.dspls)

### plotting the observed values VS predicted values for 6 components
plot(y,mod.dspls$fitted.values[,6], xlab="Observed values", ylab="Predicted values",
 main="Observed VS Predicted for 6 components")
points(-1000:1000,-1000:1000,type='l')

### plotting the regression coefficients
par(mfrow=c(3,1))

i=6
plot(1:dim(X)[2],mod.dspls$Bhat[,i],type='l',
    main=paste(" PLS , ncp =", i),
    ylab='',xlab='' )
par(oldpar)

}
\references{
H. Wold. Path Models with Latent Variables: The NIPALS Approach. In H.M. Blalock et al., editor, Quantitative Sociology: International Perspectives on Mathematical and Statistical Model Building, pages 307–357. Academic Press, 1975.
}
\author{
Louna Alsouki François Wahl
}
