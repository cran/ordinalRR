\name{summary.ordinalRR}
\alias{summary.ordinalRR}
\title{Summarize an object of class ordinalRR.}
\description{This function displays: the original call to ordinalRR(), a summary of the ordinal repeatability and reproducibility (R&R) data and its Bayesian posterior sample from JAGS, a table with point estimates of repeatability and model parameters for each rater, and a table with point estimates of R&R for each pair of raters. The estimates for repeatability and R&R are the basic/simple statistics defined in Equations (24) and (25) of Culp, Ryan, Chen, and Hamada (2018). Estimated repeatability is the proportion of matches on all pairs of repetitions for a given rater on a given part across all parts. Estimated R&R is the proportion of matches on all pairs of repetitions for a given pair of raters on a given part across all parts. On the other hand, reported parameter estimates for the De Mast-Van Wieringen Model (i.e., discrimination parameter a and cutpoints d) are the posterior medians based on a Bayesian fixed- or random-effects extension from Culp, Ryan, Chen, and Hamada (2018).}
\usage{\method{summary}{ordinalRR}(object,decimals=1,\dots)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\sQuote{ordinalRR} output from function ordinalRR().}
  \item{decimals}{\sQuote{positive integer} number of decimals to report model parameter estimates (repeatability and R&R are reported to two additional places).}
  \item{\dots}{mop up additional inputs.}
}
\author{Ken Ryan}

\references{
Culp, S.L., Ryan, K.J., Chen, J., and Hamada, M.S. (2018). \dQuote{Analysis of Repeatability and Reproducibility Studies with Ordinal Measurements.} \emph{Technometrics}, doi:10.1080/00401706.2018.1429317.

	de Mast, J. and van Wieringen, W.N. (2010). \dQuote{Modeling and Evaluating Repeatability and Reproducibility of Ordinal Classifications.} \emph{Technometrics}, 52(1), 94-106.	

}

\seealso{
\code{\link{ordinalRR}}
}
\examples{
\donttest{  
data(followup)
followup
x=preprocess(followup)
g.random<-ordinalRR(x)
summary(g.random)
}
}