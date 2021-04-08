\name{plot.features}
\alias{plot.features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot EIC
}
\description{
Plots the EICs of either level 1, 2, 3, target, or aligned features (multi-sample analysis).
}
\usage{
plot.features(dir, featureTable, data, plot.type, plotmz.tol = 0.01, plotrt.tol = 60, smooth = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
Directory to save the plots to.
}
  \item{featureTable}{
Feature table dataframe generated in previous steps containing level information.
}
  \item{data}{
"MSdata" object generated in previous steps.
}
  \item{plot.type}{
The type of features to plot. 1 for level 1 features, 2 for level 2 features, 3, for level 3 features, "target" for target features, and "aligned" for aligned features in multi-sample analysis.
}
  \item{plotmz.tol}{
Mass range in dalton for plotting EIC. Default = 0.01.
}
  \item{plotrt.tol}{
Retention time range in seconds for plotting EIC. Default = 60.
}
  \item{smooth}{
Smoothing level for plotting EIC. Default = 2.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns nothing, only plots all EICs of specified level to local disk.
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

plotdir <- "X:/Users/Sam_Shen/IPAtest_20210330"
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 1, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 2, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 3, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = "target", smooth = 2)
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
