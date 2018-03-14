\name{hist.ordinalRR}
\alias{hist.ordinalRR}
\title{Histogram for the latent part distributions from a Bayesian ordinal R&R analysis.}
\description{
This accepts an input of class \sQuote{ordinalRR} and constructs a plot of overlaid density estimates for the latent part qualities X[i] for i=1,2,...,I,I+1. Part I+1 is a new/hypothetical part with the standard normal distribution. For an example, run this package's demo or see Figure 4 from Culp, Ryan, Chen, and Hamada (2018).
}
\usage{
\method{hist}{ordinalRR}(x,x.low=-4,x.high=4,col="grey",...)
}

\arguments{
  \item{x}{\sQuote{ordinalRR} output from function ordinalRR().}
  \item{x.low}{\sQuote{scalar} lower bound for density support.}
  \item{x.high}{\sQuote{scalar} upper bound for density support.}
  \item{col}{\sQuote{vector} colors for density curves (can be a vector of color names or numbers of length 1 or I+1).}
  \item{\dots}{mop up additional inputs.}
}
\author{
Ken Ryan
}
\seealso{
  \code{\link{ordinalRR}}
}
\references{
Culp, S.L., Ryan, K.J., Chen, J., and Hamada, M.S. (2018). \dQuote{Analysis of Repeatability and Reproducibility Studies with Ordinal Measurements}. \emph{Technometrics}, doi:10.1080/00401706.2018.1429317.
}

\examples{
   \donttest{
      demo(ordinalRR)
   }
}