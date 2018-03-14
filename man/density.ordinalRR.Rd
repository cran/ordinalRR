\name{density.ordinalRR}
\alias{density.ordinalRR}
\title{Plot densities of the latent part distributions.}
\description{
This plots density estimates of the posterior distributions for repeatabilities, R&Rs, and proportions due to repeatability. This corresponds to Figure 5 from Culp, Ryan, Chen, and Hamada (2018). There will be a dark posterior predictive curve only if a random model was fit.
}
\usage{
\method{density}{ordinalRR}(x,plt.type=c("repeat", "rr", "prop", "all"),m=0, \dots)
}
\arguments{
  \item{x}{\sQuote{ordinalRR} output from function ordinalRR().}
  \item{plt.type}{\sQuote{character} plot type (must be \sQuote{repeat}, \sQuote{rr}, \sQuote{prop}, or \sQuote{all}). }
  \item{m}{\sQuote{nonnegative integer} relaxes the match criterion to within plus or minus m (see Equations (26) and (27) from Culp, Ryan, Chen, and Hamada, 2018).}
  \item{\dots}{mop up additional inputs.}
}
\author{
Ken Ryan
}
\seealso{
\code{\link{ordinalRR}}
}


\references{
Culp, S.L., Ryan, K.J., Chen, J., and Hamada, M.S. (2018). \dQuote{Analysis of Repeatability and Reproducibility Studies with Ordinal Measurements.} \emph{Technometrics}, doi:10.1080/00401706.2018.1429317.
}