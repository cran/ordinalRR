\name{preprocess}
\alias{preprocess}

\title{Format an ordinal R&R data frame into object required by function ordinalRR.}
\description{This function accepts a concise data frame containing the ordinal responses from an R&R study and expands it into the necessary multinomial counts, so it can be an input to function ordinalRR() for Bayesian data analysis.}

\usage{
preprocess(x, J=3, K=2, H=4)
}
\arguments{
  \item{x}{\sQuote{data.frame} containing J*K columns and entries from the H-point ordinal scale 1:H. The required format is a row for each part and blocks of K adjacent columns for the repetitions of each of J raters, e.g., rater 1's columns are the first K and rater J's columns are the last K.}
  \item{J}{\sQuote{positive integer} defining the number of raters.}
  \item{K}{\sQuote{positive integer} defining the number of repetitions per rater.}
  \item{H}{\sQuote{positive integer} defining the H-point ordinal scale on \{1,...,H\}.}
}
\value{
  \item{I}{\sQuote{positive integer} number of parts.}
  \item{J}{\sQuote{positive integer} number of raters.}
  \item{K}{\sQuote{positive integer} number of repetitions per rater.}
  \item{H}{\sQuote{positive integer} length of the H-point ordinal scale on \{1,...,H\}.}
  \item{x}{\sQuote{data.frame} containing J*K columns and entries from the H-point ordinal scale 1:H. Each part is a row, and there are blocks of K adjacent columns for the repetitions of each of J raters, e.g., rater 1's columns are the first K and rater J's columns are the last K.}
  \item{R}{\sQuote{array} x is expanded into a 3-dimensional array (i.e., part 1:I, operator 1:J, ordinal value 1:H) with multinomial counts.}
 \item{preprocess}{\sQuote{Boolean} will be TRUE if the data are ready for input into function ordinaRR() for Bayesian analysis with JAGS.}}

\author{Ken Ryan}

\seealso{\code{\link{followup}}}
\examples{
data(followup)
followup
preprocess(followup) 
}