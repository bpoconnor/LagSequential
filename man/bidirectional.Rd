\name{bidirectional}
\alias{bidirectional}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{bidirectional}
\description{Tests for bidirectional dependence between pairs of sequential transitions.}
\usage{
bidirectional(data, ncodes = 6, labels = NULL, lag = 1,
  adjacent = 1, tailed = 1, permtest = 0, nperms = 10,
  nblocks = 3, confid = 95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  \code{}A one-column, all-numeric dataframe or vector of code sequences, or a square
  frequency transition matrix. If data is not a frequency transition matrix,
  then data must be a series of integer codes with values ranging from "1" to what 
  ever value the user specifies in the "ncodes" argument. Cases with missing values 
  are not permitted in the data file.
}
  \item{ncodes}{
  \code{}The number of possible code values.
  }
  \item{labels}{
  \code{}Optional argument for providing labels to the code values. Accepts 
  a list of string variables. Labels should not exceed five characters. If 
  unspecified, codes will be labeled "Code1", "Code2", etc.
  }
  \item{lag}{
  \code{}The lag number for the analyses.
  }
  \item{adjacent}{
  \code{}Can adjacent values be coded the same? Enter "0" if adjacent events 
  can never be the same. Enter "1" if adjacent events can always be the same. 
  Enter "2" if some adjacent events can, and others cannot, be the same; then 
  enter the appropriate onezero matrix for your data.
  }
  \item{tailed}{
  \code{}Specify whether significance tests are one-tailed or two-tailed. 
  Options are "1" or "2".
  }
  \item{permtest}{
  \code{}Do you want to run permutation tests of significance? Options are 
  "0" for no, or "1" for yes. Warning: these computations can be time consuming.
  }
  \item{nperms}{
  \code{}The number of permutations per block.
  }
  \item{nblocks}{
  \code{}The number of blocks of permutations.
  }
  \item{confid}{
  \code{}The percentage confidence intervals for the permutation tests. Options 
  are "95" for 95 percent confidence intervals or 99 for 99 percent confidence intervals.
  }
}
\details{
This function tests the bidirectional dependence of behaviours \emph{i} to \emph{j}, 
and \emph{j} to \emph{i}, an additive sequential pattern described by Wampold 
and Margolin (1982) and Wampold (1989, 1992). Bidirectional dependence suggests 
a reciprocal effect of behaviours. That is, behaviour \emph{i} influences 
behaviour \emph{j} and behaviour \emph{j} influences behaviour \emph{i}. For 
example, if behaviour \emph{i} is a husband's positive behaviour, and behaviour 
\emph{j} is his wife's positive behaviour, a test of bidirectional dependence 
asks whether the husband reciprocates the wife's positive behaviour, \emph{and} 
the wife reciprocates the husband's positive behaviour (See Margolin and Wampold, 1982). 
Bidirectional dependence is sometimes called a "circuit".
}

\value{ A list with the following elements:

  \item{freqs}{The transitional frequency matrix}

  \item{bifreqs}{The bidirectional frequencies}

  \item{expbifreqs}{The expected bidirectional frequencies}

  \item{kappas}{The bidirectional kappas}

  \item{z}{The z values for the kappas}

  \item{pk}{The p values (significance levels) for the kappas}
}

\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for analyzing 
 	lag-sequential categorical data. \emph{Behavior Research Methods, 
 		Instrumentation, and Computers, 31,} 718-726.}
\cr\cr {Wampold, B. E., & Margolin, G. (1982). Nonparametric strategies to test 
	the independence of behavioral states in sequential data. \emph{Psychological 
		Bulletin, 92,} 755-765.}
\cr\cr {Wampold, B. E. (1989). Kappa as a measure of pattern in sequential data. 
	\emph{Quality & Quantity, 23,} 171-187.}
\cr\cr {Wampold, B. E. (1992). The intensive examination of social interactions. 
	In T. Kratochwill & J. Levin (Eds.), \emph{Single-case research design and 
		analysis: New directions for psychology and education} (pp. 93-131). 
		Hillsdale, NJ: Erlbaum.}
}
\author{Zakary A. Draper & Brian P. O'Connor}
%\note{  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{  ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
bidirectional(data_Wampold_1982, ncodes = 6,
              labels=c('HPos','HNeu','HNeg','WPos','WNeu','WNeg') )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sequential Analysis }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
