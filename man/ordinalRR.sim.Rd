\name{ordinalRR.sim}
\alias{ordinalRR.sim}

\title{Simulate an ordinal R&R data set.}
\description{This function accepts design parameters and model parameters for the random-effects model from Section 3 of Culp, Ryan, Chen, and Hamada (2018) and outputs a simulated R&R data set in preprocessed form.}

\usage{
ordinalRR.sim(H=4L,I=30L,J=3L,K=2L,mu.a=2.6,sigma.a=.2,lambda=c(11,44,29,40),seed=10L)
}
\arguments{
  \item{H}{\sQuote{positive integer} length of the H-point ordinal scale on \{1,...,H\}.}
  \item{I}{\sQuote{positive integer} number of parts.}
  \item{J}{\sQuote{positive integer} number of raters.}
  \item{K}{\sQuote{positive integer} number of repetitions per rater.}
  \item{mu.a}{\sQuote{scalar} mean of log(alpha_j).}
  \item{sigma.a}{\sQuote{positive scalars} standard deviation of log(alpha_j).}
  \item{lambda}{\sQuote{vector of length H of positive scalars} Dirichlet parameters used to induce a distribution on cutpoints.}
  \item{seed}{\sQuote{positive integer} used to set R's seed for pseudo-random number generation.}

}
\value{
  \item{I}{\sQuote{positive integer} number of parts.}
  \item{J}{\sQuote{positive integer} number of raters.}
  \item{K}{\sQuote{positive integer} number of repetitions per rater.}
  \item{H}{\sQuote{positive integer} length of the H-point ordinal scale on \{1,...,H\}.}
  \item{x}{\sQuote{data.frame} containing J*K columns and entries from the H-point ordinal scale 1:H. Each part is a row, and there are blocks of K adjacent columns for the repetitions of each of J raters, e.g., rater 1's columns are the first K and rater J's columns are the last K.}
  \item{R}{\sQuote{array} x is expanded into a 3-dimensional array (i.e., part 1:I, operator 1:J, ordinal value 1:H) with multinomial counts.}
 \item{preprocess}{\sQuote{Boolean} will be TRUE if the data are ready for input into function ordinaRR() for Bayesian analysis with JAGS.}
}

\author{Ken Ryan}

\references{
Culp, S.L., Ryan, K.J., Chen, J., and Hamada, M.S. (2018). \dQuote{Analysis of Repeatability and Reproducibility Studies with Ordinal Measurements.} \emph{Technometrics}, doi:10.1080/00401706.2018.1429317.
}

\seealso{\code{\link{preprocess}}}
\examples{
ordinalRR.sim()

}
