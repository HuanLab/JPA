\name{adduct.isotope.annotation}
\alias{adduct.isotope.annotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CAMERA Annotation
}
\description{
Performs CAMERA adduct and isotope annotation analysis.
}
\usage{
adduct.isotope.annotation(featureTable, data, polarity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{featureTable}{
Feature table dataframe generated in previous steps.
}
  \item{data}{
"MSdata" object generated in previous steps.
}
  \item{polarity}{
Polarity of the analysis, either "positive" or "negative"
}
}
\details{
CAMERA annotation can only be used for single-sample analysis.
}
\value{
Returns the input feature table with additional columns containing adduct and isotope annotation.
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
featureTable <- adduct.isotope.annotation(featureTable = featureTable, data = MSdata, polarity = "negative")
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
