\name{seqgroups}
\alias{seqgroups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{seqgroups}
\description{Computes a variety of basic sequential analysis statistics for data that are in multiple groups, or segments.}
\usage{
seqgroups(alldata, ncodes = 6, labels = NULL, lag = 1, adjacent = 1,
  onezero = NULL, tailed = 2, test = 1, output = 2, outfile = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alldata}{
  \code{}An all-numeric dataframe or vector. The data are assumed to be a series of integer codes with values ranging from "1" to what ever value the user specifies in the "ncodes" argument. The beginning of each segment, or group should be specified with a number greater than 999. Cases with missing values are not permitted in the data file.
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
  \code{}Can adjacent values be coded the same? Enter "0" if adjacent events can never be the same. Enter "1" if adjacent events can always be the same. Enter "2" if some adjacent events can, and others cannot, be the same; then enter the appropriate onezero matrix for your data.
  }
  \item{onezero}{
  \code{}Optional argument for specifying the one-zero matrix for the data. Accepts a square matrix of ones and zeros with length ncodes. A "1" indicates that the exprected frequency for a given cell is to be estimated, whereas a "0" indicates that the expected frequency for the cell should NOT be estimated, typically because it is a structural zero (codes that cannot follow one another). By default, the matrix that is created by the above commands has zeros along the main diagonal, and ones everywhere else, which will be appropriate for most data sets. However, if yor data happen to involve structural zeros that occur in cells other than the cells along the main diagonal, then you must create a ONEZERO matrix with ones and zeros that is appropriate for your data.
  }
  \item{tailed}{
  \code{}Specify whether significance tests are one-tailed or two-tailed. Options are "1" or "2".
  }
  \item{test}{
  \code{}Specify whether to run tests for homogeneity of homogeneity or stationarity. Homogeneity should be tested when groups in the data are actually different groups, whereas stationarity should be tested when groups in the data are segments of a single stream of observations. Options are "1" for homogeneity of "2" for stationarity.
  }
  \item{output}{
  \code{}Specify the desired output. Options are "1" for pooled data only, or "2" for all data sets.
  }
  \item{outfile}{
  \code{}Specify which results should be saved in an outfile.
  }
}
\details{
Simple function for computing a variety of sequential statistics for data that are in multiple groups, or segments. This is the same as the sequential() provided in this package, but allows for the data to be segmented. Sequential statistics are calculated for each segment, as well as for the data pooled across all segments.
}
\value{
For each of the groups or segments and for the pooled data, returns square matrices (the cells of which correspond with the cells of the transitional frequency matrix) for each of the following: transitional frequencies with row and column totals, expected frequencies, transitional probabilities, adjusted residuals and significance levels, Yule's Q values, transformed Kappas (Wampold, 1989, 1992, 1995), z values for the kappas, and significance levels.

  \item{Transitional frequencies }{Each cell indicates the number of times a state is followed by another state. The cell in row \emph{i} and column \emph{j} indicates the number of times state \emph{j} is followed by state \emph{i}.}

  \item{Expected frequencies }{Each cell indicates the number of times a state would be expected to follow another state based on chance alone. For example, the cell in row \emph{i} and column \emph{j} indicates the number of times state \emph{j} would be expected to follow state \emph{i} based on chance.}

  \item{Transitional probabilities }{Each cell indicates the probability that the state indicated by the column follows the state indicated by the row.}

  \item{Chi-square test }{Significance test of difference between the transitional frequencies matrix and the expected frequencies matrix.}

  \item{Adjusted Residuals }{Given by the observed count minus the expected count divided by an estimate of the standard error. Adjusted residuals are standardized values with a mean of 0 and standard deviation of 1. }

  \item{Significance sevels }{Significance of difference between transitional and expected frequencies for each cell. }

  \item{Yule's Q values }{An indication of the strength of the relationship between the antecedent and the consequence (i.e., an index of effect size). }

  \item{Transformed kappas }{Kappa provides an index of effect size. Values range from -1 to 1, with 0 being no effect. These values may be used in subsequent analyses, as they represent the extent to which the antecedent predicts the consequence for each case. }

  \item{z values for kappas }{Standardized kappa values. }

  \item{Significance levels of kappas }{p-values for kappas. }
}
\references{
 {O'Connor, B. P. (1999). Simple and flexible SAS and SPSS programs for analyzing lag-sequential categorical data. Behavior Research Methods, Instruments & Computers, 31, 718-726.}
\cr\cr {Wampold, B. E. (1989). Kappa as a measure of pattern in sequential data. Quality & Quantity, 23, 171-187.}
\cr\cr {Wampold, B. E. (1992). The intensive examination of social interactions. In T. Kratochwill & J. Levin (Eds.), Single-case research design and analysis: New directions for psychology and education (pp. 93-131). Hillsdale, NJ: Erlbaum.}
\cr\cr {Wampold, B. E. (1995). Analysis of behavior sequences in psychotherapy. In J. Siegfried (Ed.), Therapeutic and everyday discourse as behavior change: Towards a micro-analysis in psychotherapy process research (pp. 189-214). Norwood, NJ: Ablex.}
}
\author{Zakary A. Draper & Brian P. O'Connor}
%\note{  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{  ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
seqgroups(data_seqgroups, ncodes = 6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sequential Analysis }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
