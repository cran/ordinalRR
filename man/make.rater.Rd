\name{make.rater}
\alias{make.rater}

\title{Format the parameters for a single rater.}
\description{This function accepts fixed parameters for the De Mast-Van Wieringen Model of a single rater and formats them, so meaningful parametric functions can be displayed with the plot() function.}

\usage{
make.rater(alpha,cutpoints)
}
\arguments{
  \item{alpha}{\sQuote{positive scalar} the rater's discrimination parameter.}
  \item{cutpoints}{\sQuote{vector} H-1 (ordered values) corresponding to an H-point ordinal scale.}
}
\value{
\item{rater}{\sQuote{rater} a list of rater parameters}

}

\author{Ken Ryan}

\seealso{\code{\link{compute.q}}}

\references{
	de Mast, J. and van Wieringen, W.N. (2010). \dQuote{Modeling and Evaluating Repeatability and Reproducibility of Ordinal Classifications.} \emph{Technometrics}, 52(1), 94-106.	
}


\examples{
(rater1=make.rater(1,c(-1.7,-0.5,1.6))) #3 cutpoints, so H=4
(prob=compute.q(rater1,.1)) #probabilities of 1,2,3,4 if part quality is x=.1
sum(prob) #should sum to one
}


