\name{twocells}
\alias{twocells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{twocells}
\description{Simultaneously tests the the unidirectional dependence of \emph{i} to \emph{j}, and the unidirectional dependence of \emph{k} to \emph{L}, an additive pattern described by Wampold and Margolin (1982) and Wampold (1989, 1992).}
\usage{
twocells(data, i, j, k, L, ncodes = 6, labels = NULL, lag = 1,
  adjacent = 1, tailed = 1, permtest = 0, nperms = 10,
  nblocks = 3, confid = 95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe or vector. The data are assumed to be a series of integer codes with values ranging from "1" to what ever value the user specifies in the "ncodes" argument.

         Cases with missing values are not permitted in the data file.
  }
  \item{i}{
  \code{}Code value for \emph{i}.
  }
  \item{j}{
  \code{}Code value for \emph{j}.
  }
  \item{k}{
  \code{}Code value for \emph{k}.
  }
  \item{L}{
  \code{}Code value for \emph{L}.
  }
  \item{ncodes}{
  \code{}The number of possible code values.
  }
  \item{labels}{
  \code{}Optional argument for providing labels to the code values. Accepts a list of string variables. Labels should not exceed five characters. If unspecified, codes will be labeled "Code1", "Code2", etc.
  }
  \item{lag}{
  \code{}The lag number for the analyses.
  }
  \item{adjacent}{
  \code{}Can adjacent values be coded the same? Options are "1" for yes, and "0" for no.
  }
  \item{tailed}{
  \code{}Specify whether significance tests are one-tailed or two-tailed. Options are "1" or "2".
  }
  \item{permtest}{
  \code{}Do you want to run permutation tests of significance? Options are "0" for no, or "1" for yes. Warning: these computations can be time consuming.
  }
  \item{nperms}{
  \code{}The number of permutations per block.
  }
  \item{nblocks}{
  \code{}The number of blocks of permutations.
  }
  \item{confid}{
  \code{}The percentage confidence intervals for the permutation tests. Options are "95" for 95 percent confidence intervals or 99 for 99 percent confidence intervals.
  }
}
\details{
This function simultaneously tests the unidirectional dependence of \emph{i} to \emph{j} and the unidirectional dependence of \emph{k} to \emph{L}. The user specifies the code values used for \emph{i}, \emph{j}, \emph{k}, and \emph{L} in the analyses. This function is useful for answering theoretical questions about behaviour. For example, Wampold and Margolin (1982) described a situation wherein a spouse responds to negative behaviours with something other than a negative behaviour.
}
\value{
Returns the transitional frequency matrix, observed and expected values for the two cell test, kappa, the z value for kappa, and the significance level.

  \item{Transitional frequncies}{Each cell indicates the number of times a state is followed by another state. The cell in row \emph{i} and column \emph{j} indicates the number of times state \emph{j} is followed by state \emph{i}.}

  \item{Expected frequency}{The expected number of transitions from \emph{i} to \emph{j} and from \emph{k} to \emph{L} based on chance alone.}

  \item{Observed frequency}{The observed number of transitions from \emph{i} to \emph{j} and from \emph{k} to \emph{L}.}

  \item{Kappa}{Kappa provides an index of effect size. Values range from -1 to 1, with 0 being no effect. These values may be used in subsequent analyses.}

  \item{z value for kappa}{Standardized kappa value.}

  \item{Significance level of kappa}{p-values for kappa.}

}
\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for analyzing lag-sequential categorical data. Behavior Research Methods, Instruments & Computers, 31, 718-726.}
\cr\cr {Wampold, B. E. (1989). Kappa as a measure of pattern in sequential data. Quality & Quantity, 23, 171-187.}
\cr\cr {Wampold, B. E. (1992). The intensive examination of social interactions. In T. Kratochwill & J. Levin (Eds.), Single-case research design and analysis: New directions for psychology and education (pp. 93-131). Hillsdale, NJ: Erlbaum.}
\cr\cr {Wampold, B. E., & Margolin, G. (1982). Nonparametric strategies to test the independence of behavioral states in sequential data. Psychological Bulletin, 92, 755-765.}
}
\author{Zakary A. Draper & Brian P. O'Connor}
%\note{  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{  ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
twocells(data_sequential, ncodes = 6, i = 1, j = 6, k = 5, L = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sequential Analysis }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
