\name{XCMS.featureTable}
\alias{XCMS.featureTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate XCMS Feature Table
}
\description{
Generate feature table from raw .mzXML files using XCMS's CentWave peak detection algorithm.
}
\usage{
XCMS.featureTable <- function(dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01, snthresh = 6, integrate = 1,
                         prefilter = c(3,100), noise = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
Directory of folder containing all .mzXML files of interest. This folder must contain only relevant .mzXML files of DDA sample(s).
}
  \item{mz.tol}{
Mass tolerance in ppm used to identify level 1 features by matching the feature's mz against all ms2 precursor mz. Default = 10.
}
  \item{ppm}{
Numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition. Default = 10.
}
  \item{peakwidth}{
Numeric(2) with the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds. Default = c(5,20).
}
  \item{mzdiff}{
Numeric(1) representing the minimum difference in m/z dimension required for peaks with overlapping retention times; can be negative to allow overlap. During peak post-processing, peaks defined to be overlapping are reduced to the one peak with the largest signal. Default = 0.01.
}
  \item{snthresh}{
Numeric(1) defining the signal to noise ratio cutoff. Default = 6.
}
  \item{integrate}{
Integration method. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact. Default = 1.
}
  \item{prefilter}{
Numeric(2): c(k, I) specifying the prefilter step for the first analysis step (ROI detection). Mass traces are only retained if they contain at least k peaks with intensity >= I. Default = c(3,100).
}
  \item{noise}{
Numeric(1) allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection). Default = 100.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a feature table in dataframe format, containing the following columns in order: "mz", "rt", "rtmin", "rtmax", "maxo", "sample", "level" where "maxo" is the feature peak intensity. Also returns an XCMS object that is used by functions later in the IPA workflow.
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
