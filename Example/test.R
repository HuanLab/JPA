library(IPA)
# Example workflow scripts
#########################################################################################
# Workflow 1: XCMS -> level3 -> alignment -> plot -> annotation
# Single
dir = "X:/Users/Sam_Shen/IPAtest_20210330/singleDDA"
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
featureTable <- find.level3features(data = MSdata)

# Multi
dir = "X:/Users/Sam_Shen/IPAtest_20210330/multiDDA"
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
featureTable <- find.level3features(data = MSdata)
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
#########################################################################################
# Workflow 2: XCMS -> level3 -> target -> alignment -> plot -> annotation
# Single
dir = "X:/Users/Sam_Shen/IPAtest_20210330/singleDDA"
tarFTdir = "X:/Users/Sam_Shen/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
featureTable <- find.level3features(data = MSdata)
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)

# Multi
dir = "X:/Users/Sam_Shen/IPAtest_20210330/multiDDA"
tarFTdir = "X:/Users/Sam_Shen/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
featureTable <- find.level3features(data = MSdata)
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
#########################################################################################
# Workflow 3: Existing Feature table -> level3 -> alignment -> plot -> annotation
# Single
dir = "X:/Users/Sam_Shen/IPAtest_20210330/singleDDA"
FTdir = "X:/Users/Sam_Shen/IPAtest_20210330/singleDDA"
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
featureTable <- find.level3features(data = MSdata)

# Multi
dir = "X:/Users/Sam_Shen/IPAtest_20210330/multiDDA"
FTdir = "X:/Users/Sam_Shen/IPAtest_20210330/multiDDA"
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
featureTable <- find.level3features(data = MSdata)
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
#########################################################################################
# Workflow 4: Existing Feature table -> level3 -> target -> alignment -> plot -> annotation
# Single
dir = "X:/Users/Sam_Shen/IPAtest_20210330/singleDDA"
FTdir = "X:/Users/Sam_Shen/IPAtest_20210330/singleDDA"
tarFTdir = "X:/Users/Sam_Shen/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
featureTable <- find.level3features(data = MSdata)
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)

# Multi
dir = "X:/Users/Sam_Shen/IPAtest_20210330/multiDDA"
FTdir = "X:/Users/Sam_Shen/IPAtest_20210330/multiDDA"
tarFTdir = "X:/Users/Sam_Shen/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
featureTable <- find.level3features(data = MSdata)
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
#########################################################################################

# Continuation of Workflows 1, 2, 3, 4
plotdir <- "X:/Users/Sam_Shen/IPAtest_20210330"
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 1, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 2, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 3, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = "target", smooth = 2)
# only for multi-sample analysis
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = "multi", smooth = 2)

lib_directory <- "E:/SAM"
lib_name <- "convertedLibraryNeg.msp"
featureTable <- ms2.tofeaturetable(data = MSdata, featureTable = featureTable)
featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.1)
# CAMERA adduct and isotope annotation only available for single sample analysis
featureTable <- adduct.isotope.annotation(featureTable = featureTable, data = MSdata, polarity = "negative")





