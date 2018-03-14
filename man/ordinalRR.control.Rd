\name{ordinalRR.control}
\alias{ordinalRR.control}

\title{Set control parameters for a Bayesian ordinal R&R model.}
\description{
The R-function \sQuote{ordinalRR.control} sets various control parameters for the prior if using the random-effects version of the model and for the MCMC. See Section 3.2 of Culp, Ryan, Chen, and Hamada (2018) for more on this prior. Default settings match this paper.
}
\usage{
ordinalRR.control(mu.mu.alpha = 0.8, tau.mu.alpha = 0.4, mu.tau.alpha = 4,
tau.tau.alpha = 0.4, mu.lambda = 2, tau.lambda = 0.2, rjags.B = 10000L,
rjags.Burn = 1000L, rjags.n.chains = 1L, rjags.n.adapt = 5000L,r.seed=10L,rjags.seed=10L)
}
\arguments{
  \item{mu.mu.alpha}{\sQuote{positive scalar} mean of the normal prior for mu.alpha.}
  \item{tau.mu.alpha}{\sQuote{positive scalar} precision=1/variance of the normal prior for mu.alpha.}
  \item{mu.tau.alpha}{\sQuote{positive scalar} mean of the log-normal prior for tau.alpha.}
  \item{tau.tau.alpha}{\sQuote{positive scalar} precision of the log-normal prior for tau.alpha.}
  \item{mu.lambda}{\sQuote{positive scalar} mean of the log-normal prior for the lambda.h.}
  \item{tau.lambda}{\sQuote{positive scalar} precision of the log-normal prior for the lambda.h.}
  \item{rjags.B}{\sQuote{positive integer} length of JAGS MCMC chain retained.}
  \item{rjags.Burn}{\sQuote{positive integer} length of initial JAGS MCMC chain burnin discarded.}
  \item{rjags.n.chains}{\sQuote{1} number of JAGS MCMC chains (currently only programmed to accept 1).}
  \item{rjags.n.adapt}{\sQuote{positive integer} rjags n.adapt parameter within command jag.model().}
  \item{r.seed}{\sQuote{positive integer} sets seed within R during function ordinalRR(). This is for predictive inference on a new rater which is only used with the random-effects model. This does not fix the JAGS seed for the MCMC.}
  \item{rjags.seed}{\sQuote{positive integer} sets seed within JAGS for the posterior sample from function ordinalRR().}
  

}

\author{Ken Ryan}

\references{
Culp, S.L., Ryan, K.J., Chen, J., and Hamada, M.S. (2018). \dQuote{Analysis of Repeatability and Reproducibility Studies with Ordinal Measurements.} \emph{Technometrics}, doi:10.1080/00401706.2018.1429317.

Plummer, M. (2016). \dQuote{RJAGS: Bayesian Graphical Models using MCMC.} R Package Version 4-6, \url{https://CRAN.R-project.org/package=rjags}. 

Plummer, M. (2017). \dQuote{JAGS: A Program for Analysis of Bayesian Graphical Models using Gibbs Sampling.} Version 4.3.0, \url{http://mcmc-jags.sourceforge.net}.
}

\keyword{classes}
\keyword{models}
\keyword{methods}