\name{ms2.tofeaturetable}
\alias{ms2.tofeaturetable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assigns MS2 to Corresponding Features
}
\description{
Extracts and assigns MS2 spectra to features in feature table.
}
\usage{
ms2.tofeaturetable(data, featureTable, rt.tol = 30, mz.tol = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
"MSdata" object generated in previous steps.
}
  \item{featureTable}{
Feature table dataframe generated in previous steps.
}
  \item{rt.tol}{
Retention time tolerance (seconds) used in MS2 assignment to feature table. Default = 30.
}
  \item{mz.tol}{
Mass tolerance tolerance (Da) used in MS2 assignment to feature table. Default = 0.01.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a feature table with additional columns containing MS2 fragment m/z and intensity information for each feature.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sam Shen, Jian Guo, Tao Huan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(IPA)
dir = "X:/Users/Sam_Shen/IPAtest_20210330/singleDDA"
featureTable <- peak.picking(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
featureTable <- find.level3features(data = MSdata)
featureTable <- ms2.tofeaturetable(data = MSdata, featureTable = featureTable)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
