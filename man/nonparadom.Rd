\name{nonparadom}
\alias{nonparadom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{nonparadom}
\description{Tests for nonparallel dominance, a form of asymmetry in predictability, between i to j and k to L (Wampold, 1984, 1989, 1992, 1995).}
\usage{
nonparadom(data, ncodes = 6, i, j, k, L, labels = NULL, lag = 1, adjacent = 1, tailed = 1, permtest = 0, nperms = 10, nblocks = 3, confid = 95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe or vector. The data are assumed to be a series of integer codes with values ranging from "1" to what ever value the user specifies in the "ncodes" argument. Cases with missing values are not permitted in the data file.
}
  \item{ncodes}{
  \code{}The number of possible code values.
  }
  \item{i}{
  \code{}Code value for i.
  }
  \item{j}{
  \code{}Code value for j.
  }
  \item{k}{
  \code{}Code value for k.
  }
  \item{L}{
  \code{}Code value for L.
  }
  \item{labels}{
  \code{}Optional argument for providing labels to the code values. Accepts a list of string variables. Labels should not exceed five characters. If unspecified, codes will be labeled "Code1", "Code2", etc.
  }
  \item{lag}{
  \code{}The lag number for the analyses.
  }
  \item{adjacent}{
  \code{}Can adjacent values be coded the same? Options are "1" for yes or "0" for no.
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
Tests for nonparallel dominance or asymmetry in predictability, which is the difference in predictability between i to j and k to L (e.g., whether a partner follows a negative behaviour with a positive behaviour), as described by Wampold (1984, 1989, 1992, 1995).
}
\value{
Returns the transitional frequency matrix, expected fre- quencies, expected and observed nonparallel dominance frequencies, kappas, the z value for kappas, and significance levels.

  \item{Transitional frequencies}{Each cell indicates the number of times a state is followed by another state. The cell in row \emph{i} and column \emph{j} indicates the number of times state \emph{j} is followed by state \emph{i}.}

  \item{Expected frequencies}{Each cell indicates the number of times a state would be expected to follow another state based on chance alone. For example, the cell in row \emph{i} and column \emph{j} indicates the number of times state \emph{j} would be expected to follow state \emph{i} based on chance.}

  \item{Dominance frequencies}{}

  \item{Expected dominance frequencies}{}

  \item{Dominance types}{There are 4 sequential dominance case types described by Wampold (1989). These cases describe the direction of the effect for \emph{i} on \emph{j} and \emph{j} on \emph{i}. The four cases are: (1) \emph{i} increases \emph{j}, and \emph{j} increases \emph{i}, (2) \emph{i} decreases \emph{j}, and \emph{j} decreases \emph{i}, (3) \emph{i} increases \emph{j}, and \emph{j} decreases \emph{i}, and (4) \emph{i} decreases \emph{j}, and \emph{j} increases \emph{i}. Each cell of this matrix indicates the case that applies to the transition indicated by the cell.}

  \item{Dominance kappas}{An index of effect size for nonparallel dominance. Kappa ranges from -1 to 1, with 0 being no effect.}

  \item{z values for kappas}{Standardized kappa values.}

  \item{Significance levels of kappas}{p-values for kappas.}
}

\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for analyzing lag-sequential categorical data. Behavior Research Methods, Instruments & Computers, 31, 718-726.}
\cr\cr {Wampold, B. E. (1984). Tests of dominance in sequential categorical data. Psychological Bulletin, 96, 424-429.}
\cr\cr {Wampold, B. E. (1989). Kappa as a measure of pattern in sequential data. Quality & Quantity, 23, 171-187.}
\cr\cr {Wampold, B. E. (1992). The intensive examination of social interactions. In T. Kratochwill & J. Levin (Eds.), Single-case research design and analysis: New directions for psychology and education (pp. 93-131). Hillsdale, NJ: Erlbaum.}
\cr\cr {Wampold, B. E. (1995). Analysis of behavior sequences in psychotherapy. In J. Siegfried (Ed.), Therapeutic and everyday discourse as behavior change: Towards a micro-analysis in psychotherapy process research (pp. 189-214). Norwood, NJ: Ablex.}
}
\author{Zakary A. Draper & Brian P. O'Connor}
%\note{  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{  ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
nonparadom(data_sequential, ncodes = 6, i = 5, j = 3, k = 6, L = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sequential Analysis }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
