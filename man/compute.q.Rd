\name{compute.q}
\alias{compute.q}

\title{Compute the probabilities for a single rater at a fixed part quality.}
\description{This function accepts an object of class \sQuote{rater} and the latent (scalar) value of a hypothetical part's quality. The probabilities of the H ordinal categories are outputted based on the De Mast-Van Wieringen Model.}

\usage{
compute.q(rater,x)
}
\arguments{
  \item{rater}{\sQuote{rater} the rater's parameters. (Use make.rater() to construct this input.)}
  \item{x}{\sQuote{scalar} latent quality of a fixed part.}
}
\value{
\item{vector}{a probability vector of length H}
}

\author{Ken Ryan}

\seealso{\code{\link{make.rater}}}

\references{
	de Mast, J. and van Wieringen, W.N. (2010). \dQuote{Modeling and Evaluating Repeatability and Reproducibility of Ordinal Classifications.} \emph{Technometrics}, 52(1), 94-106.	
}


\examples{
(rater1=make.rater(1,c(-1.7,-0.5,1.6))) #3 cutpoints, so H=4
(prob=compute.q(rater1,.1)) #probabilities of 1,2,3,4 if part quality is x=.1
sum(prob) #should sum to one
}


