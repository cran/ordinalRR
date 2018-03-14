\name{print.ordinalRR}
\alias{print.ordinalRR}
\title{Print function for an object of class ordinalRR.}
\description{This standard print function displays: the function call to ordinalRR() that produced this object and descriptions of the ordinal R&R data and its Bayesian posterior sample from JAGS. }
\usage{\method{print}{ordinalRR}(x,\dots)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\sQuote{ordinalRR} output from function ordinalRR().}
  \item{\dots}{mop up additional inputs.}
}
\author{Ken Ryan}


\seealso{
\code{\link{ordinalRR}}
}
\examples{
\donttest{  
data(followup)
followup
x=preprocess(followup)
g.random<-ordinalRR(x)
print(g.random)
}
}