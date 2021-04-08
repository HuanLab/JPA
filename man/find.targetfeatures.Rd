\name{find.targetfeatures}
\alias{find.targetfeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find User Defined Target Features
}
\description{
Adds user defined target features to original feature table if they do not exist already.
}
\usage{
find.targetfeatures(data, tarFTdir, tarFTname, mass.tol = 0.01, rt.tol = 30, target.threshold = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
"MSdata" object generated in previous steps by functions "XCMS.featureTable", "custom.featureTable", or "find.level3features".
}
  \item{tarFTdir}{
Directory containing target feature table in CSV format.
}
  \item{tarFTname}{
Name of feature table in CSV format.
}
  \item{mass.tol}{
Mass tolerance in daltons used to determine if a target feature already exists. Default = 0.01.
}
  \item{rt.tol}{
Retention time tolerance in seconds used to determine if a target feature already exists. Default = 30.
}
  \item{target.threshold}{
Intensity threshold used to extract target features. The target feature's peak intensity must be higher than target.threshold * average intensity to be valid. Default = 2.
}
}
\details{
The provided target feature table must be in CSV format and contains only the following columns in order: "mz", "rt".
}
\value{
Returns an updated feature table in dataframe format that includes newly added target features. Also returns an updated XCMS object that is used by functions later in the IPA workflow.
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
tarFTdir = "X:/Users/Sam_Shen/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
featureTable <- peak.picking(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
featureTable <- find.level3features(data = MSdata)
featureTable <- add.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)
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
