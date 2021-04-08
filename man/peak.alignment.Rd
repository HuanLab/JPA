\name{peak.alignment}
\alias{peak.alignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Peak Alignment for Multi-sample Analysis
}
\description{
Performs XCMS peak alignment for all level 1, 2, 3, and target features identified in multi-sample analysis.
}
\usage{
peak.alignment(data, bw = 5, minfrac = 0.5, mzwid = 0.015, minsamp = 1, max = 100, quantitative.method = "maxo")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
"MSdata" object generated in previous steps.
}
  \item{bw}{
numeric(1) defining the bandwidth (standard deviation ot the smoothing kernel) to be used. Default = 5.
}
  \item{minfrac}{
numeric(1) defining the minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group (feature). Default = 0.5.
}
  \item{mzwid}{
numeric(1) defining the size of the overlapping slices in mz dimension. Default = 0.015.
}
  \item{minsamp}{
numeric(1) with the minimum number of samples in at least one sample group in which the peaks have to be detected to be considered a peak group (feature). Default = 1.
}
  \item{max}{
numeric(1) with the maximum number of peak groups to be identified in a single mz slice. Default = 100.
}
  \item{quantitative.method}{
Definition of feature intensity: either "maxo" for peak height or "into" for peak area.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns aligned feature table with an additional column containing intensity values for each additional sample in multi-sample analysis. Note: features' level information will be lost at this point.
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
dir = "X:/Users/Sam_Shen/IPAtest_20210330/multiDDA"
featureTable <- peak.picking(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
featureTable <- find.level3features(data = MSdata)
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
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
