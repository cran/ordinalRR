\name{ordinalRR}
\alias{ordinalRR}

\title{Fit a Bayesian ordinal R&R model using JAGS.}
\description{
This function can fit either the fixed- or random-effects model from Section 3 of Culp, Ryan, Chen, and Hamada (2018). The outputted class is of type \sQuote{ordinalRR}, and there are S3 generic functions (e.g., plot and density) for class \sQuote{ordinalRR} that make graphs like those from the referenced paper. The user can also use the posterior sample for a customized Bayesian data analysis; see the value list for details on how the posterior sample from JAGS is outputted.
}

\usage{
ordinalRR(x, random = TRUE, control = ordinalRR.control())
}
\arguments{
  \item{x}{\sQuote{list} output from function preprocess() containing design parameters and the multinomial counts.}
  \item{random}{\sQuote{Boolean} True for a random-effects model (Section 3.1 of Culp, Ryan, Chen, and Hamada, 2018); False for fixed effect (Section 3.2)}
  \item{control}{\sQuote{list} contains random-effect hyperparameters for the prior and JAGS MCMC parameters.}
}
\value{
\item{post}{\sQuote{mcmc.list} Bayesian model fit from JAGS including the posterior (i.e., a matrix with  samplecontrol$rjag.B rows and one column for each parameter).}	
\item{x}{\sQuote{mcmc} the columns X[1],...,X[I] from the posterior sample for the latent parts.}
\item{a}{\sQuote{mcmc} the columns alpha[1],...,alpha[J] from the posterior sample for the discrimination parameters.}
\item{d}{\sQuote{mcmc} the columns delta[1,1:(H-1)],...,delta[J,1:(H-1)] from the posterior sample for the cutpoints.}
\item{mu.a}{\sQuote{mcmc} column mu.alpha from the posterior sample for the mean of the normal log(alpha[J+1]) distribution (only outputed when random=TRUE).}
\item{sigma.a}{\sQuote{mcmc} column sigma.alpha from the posterior sample for the standard deviation of the normal log(alpha[J+1]) distribution (only outputed when random=TRUE).}
\item{lambda}{\sQuote{mcmc} columns lambda[1],...,lambda[H] from the posterior sample for the Dirichlet distribution used to induce a distribution on cutpoints delta[J+1,1:H] (only outputed when random=TRUE).}
\item{\dots}{remaining internal outputs used for necessary processing.}
}

\references{
Culp, S.L., Ryan, K.J., Chen, J., and Hamada, M.S. (2018). \dQuote{Analysis of Repeatability and Reproducibility Studies with Ordinal Measurements.} \emph{Technometrics}, doi:10.1080/00401706.2018.1429317.

Plummer, M. (2016). \dQuote{RJAGS: Bayesian Graphical Models using MCMC.} R Package Version 4-6, \url{https://CRAN.R-project.org/package=rjags}. 

Plummer, M. (2017). \dQuote{JAGS: A Program for Analysis of Bayesian Graphical Models using Gibbs Sampling.} Version 4.3.0, \url{http://mcmc-jags.sourceforge.net}.
}

\author{Ken Ryan}

\seealso{
\code{\link{preprocess}}
\code{\link{ordinalRR.control}}
}
\examples{
\donttest{          
data(followup)
followup
x=preprocess(followup)
g.random<-ordinalRR(x)
g.fixed<-ordinalRR(x,random=F)
}
\dontshow{
data(followup)
x=preprocess(followup[1:10,1:4],J=2,K=2,H=4)
g=ordinalRR(x,control=ordinalRR.control(rjags.B=10,rjags.Burn=10))
}
}
\keyword{classes}
\keyword{models}
\keyword{methods}