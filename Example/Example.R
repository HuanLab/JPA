library(IPA)
# Example workflow scripts
#########################################################################################
# Workflow 1: JPA-PP(centWave) -> JPA-MR -> alignment -> plot -> annotation
# Single file
dir = "X:/Users/IPAtest_20210330/singleDDA"
# Extract JPA-PP fetures 
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
# Extract JPA-MR fetures
featureTable <- find.level3features(data = MSdata, level3.threshold = 2)

# Multi files
dir = "X:/Users/IPAtest_20210330/multiDDA"
# Extract JPA-PP fetures
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
# Extract JPA-MR fetures
featureTable <- find.level3features(data = MSdata)
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
#########################################################################################
# Workflow 2: JPA-PP(centWave) -> JPA-MR -> JPA-TL -> alignment -> plot -> annotation
# Single
dir = "X:/Users/IPAtest_20210330/singleDDA"
tarFTdir = "X:/Users/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
# Extract JPA-PP fetures
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
# Extract JPA-MR fetures                             
featureTable <- find.level3features(data = MSdata)
# Extract JPA-TL fetures
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)

# Multi
dir = "X:/Users/Sam_Shen/IPAtest_20210330/multiDDA"
tarFTdir = "X:/Users/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
# Extract JPA-PP fetures
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
# Extract JPA-MR fetures
featureTable <- find.level3features(data = MSdata)
# Extract JPA-TL fetures
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)
# Deisotope the features
featureTable <- adduct.isotope.annotation(data = MSdata, polarity = "negative")
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
#########################################################################################
# Workflow 3: Existing Feature table -> JPA-MR -> alignment -> plot -> annotation
# Single
dir = "X:/Users/IPAtest_20210330/singleDDA"
FTdir = "X:/Users/IPAtest_20210330/singleDDA"
# Input JPA-PP fetures from other software
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
# Extract JPA-MR fetures
featureTable <- find.level3features(data = MSdata)

# Multi
dir = "X:/Users/IPAtest_20210330/multiDDA"
FTdir = "X:/Users/IPAtest_20210330/multiDDA"
# Input JPA-PP fetures from other software
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
# Extract JPA-MR fetures
featureTable <- find.level3features(data = MSdata)
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
#########################################################################################
# Workflow 4: Existing Feature table -> JPA-MR -> JPA-TL -> alignment -> plot -> annotation
# Single
dir = "X:/Users/IPAtest_20210330/singleDDA"
FTdir = "X:/Users/IPAtest_20210330/singleDDA"
tarFTdir = "X:/Users/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
# Input JPA-PP fetures from other software
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
# Extract JPA-MR fetures
featureTable <- find.level3features(data = MSdata)
# Extract JPA-TL fetures
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)

# Multi
dir = "X:/Users/IPAtest_20210330/multiDDA"
FTdir = "X:/Users/IPAtest_20210330/multiDDA"
tarFTdir = "X:/Users/IPAtest_20210330"
tarFTname = "LibraryCSVHILIC-.csv"
# Input JPA-PP fetures from other software
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
# Extract JPA-MR fetures
featureTable <- find.level3features(data = MSdata)
# Extract JPA-TL fetures
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
#########################################################################################

# Continuation of Workflows 1, 2, 3, 4
plotdir <- "X:/Users/IPAtest_20210330"
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 1, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 2, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 3, smooth = 2)
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = "target", smooth = 2)
# only for multi-sample analysis
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = "multi", smooth = 2)

lib_directory <- "X:/Users/Library"
lib_name <- "convertedLibraryNeg.msp"
featureTable <- ms2.tofeaturetable(data = MSdata, featureTable = featureTable)
featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.1)
# CAMERA adduct and isotope annotation only available for single sample analysis
featureTable <- adduct.isotope.annotation(featureTable = featureTable, data = MSdata, polarity = "negative")





