JPA R Package User Manual
================
Jian Guo, Sam Shen, Tao Huan
07/04/2021

-   [Part 1: Introduction and
    Installation](#part-1-introduction-and-installation)
-   [Part 2: MS1 peak picking Feature
    Extraction](#part-2-ms1-PP-feature-extraction)
    -   [2.1 XCMS Feature Extraction](#21-xcms-feature-extraction)
    -   [2.2 Custom Featuretable Input](#22-custom-featuretable-input)
-   [Part 3: MS2 recognition Feature
    Identification](#part-3-MR-feature-identification)
-   [Part 4: Add Target Features](#part-4-add-target-features)
-   [Part 5: Sample Alignment](#part-5-sample-alignment)
-   [Part 6: MS2 Annotation](#part-6-ms2-annotation)
-   [Part 7: Export EIC](#part-7-export-eic)
-   [Part 8: CAMERA Annotation](#part-8-camera-annotation)
-   [Part 9: Additional Details and
    Notes](#part-9-additional-details-and-notes)

# Part 1: Introduction and Installation

`JPA` is an R package for identifying MS2 recognition and targeted list features and performing
integrated feature extraction and annotation. The package is written in
the language R and its source code is publicly available at
<https://github.com/HuanLab/JPA.git>.

To install `JPA` package R version 4.0.0 or above is required, and we
recommend using RStudio to complete the installation and usage of
`JPA` by following the steps below:

``` r
# Install "BiocManager" package from CRAN if you do not already have it installed.
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
# Install "devtools" package from CRAN if you do not already have it installed.
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}

# Load "devtools" package.
library(devtools)

# Install "JPA" from Github using "devtools".
if (!requireNamespace("JPA", quietly = TRUE)){
  install_github("HuanLab/JPA")
}

# Load "JPA" package.
library(JPA)
```

# Part 2: MS1 Peak Picking Feature Extraction

`JPA` supports multiple ways to generate an MS1 featuretable and to
label features as "PP". Users can choose to use `XCMS` to
extract features from mzXML files (Section 2.1) or upload their own
featuretable in csv format generated by other commonly used data
processing software (Section 2.2). For the rest of the tutorial, to view
details and additional parameters of functions, type:
`help("<function name>")`. Note: for multi-sample analysis, sample
alignment is performed after MR and target features have been
identified. It will be discussed in section 5.

## 2.1 XCMS Feature Extraction

One or multiple mzXML files from DDA analyses can be analyzed at once
using XCMS to extract MS1 PP and MR features. All mzXML file(s) need
to be placed in a separate folder containing no other irrelevant mzXML
files. Additional details of XCMS algorithm is available at their
official website at: <https://rdrr.io/bioc/xcms/man/>.

``` r
# dir specifies the full directory of the folder containing mzXML file(s).
dir = "X:/Users/JPAtest_20210330/multiDDA"

# The XCMS.featureTable() function outputs a dataframe formatted featuretable as well as an MSnbase object.
# The "level" column shows the level of each feature.
featureTable <- XCMS.featureTable(dir = dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01,
                             snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100)
head(featureTable)
```

    ##          mz       rt    rtmin    rtmax  maxo sample level
    ## F1 44.99865 1446.592 1403.316 1461.054 84702      3    PP
    ## F2 44.99869 1447.319 1430.959 1457.881 85878      1    PP
    ## F3 44.99874 1680.279 1676.622 1691.061 68456      1    PP
    ## F4 44.99877 1679.435 1675.546 1691.578 66910      2    PP
    ## F5 44.99878 1680.822 1676.406 1691.724 75678      3    PP
    ## F6 44.99880 1443.377 1435.238 1455.427 83670      2    PP

## 2.2 Custom Featuretable Input

When using a custom featuretable (eg. from MS-DIAL, MZmine2, etc) for
`JPA` analyses, in order for `JPA` to successfully read the provided CSV
file, it must contain only columns in the following order: m/z,
retention time, min retention time, max retention time, intensity. For
multi-sample analysis, a CSV feature table needs to be provided for each
single mzXML file to be analyzed. Aligned feature table is not
acceptable. All CSV file(s) of interest need to be placed in a separate
folder containing no other irrelevant CSV files. Note: column 3 and
column 4 are the retention time of the feature edges, and all three
columns containing retention time information should be in seconds.

``` r
# ft_directory specifies the directory of the custom csv file(s).
FTdir = "X:/Users/JPAtest_20210330/multiDDA"
# dir specifies the full directory of the folder containing mzXML file(s).
dir = "X:/Users/JPAtest_20210330/multiDDA"

# Sample csv feature table.
setwd(FTdir)
head(read.csv(list.files(pattern = ".csv")[1], header = T, stringsAsFactors = F))
```

    ##          mz     rt     rtmin    rtmax   Height
    ## 1 194.90540 25.444 19.459002 35.47100 5737.125
    ## 2  85.02976 29.384  7.153002 70.32402 1118.750
    ## 3  87.00925 29.384  7.531002 70.32402 2168.750
    ## 4 111.00890 29.384  7.153002 72.63300 4777.000
    ## 5 173.00910 29.384  4.920000 75.61200 1909.250
    ## 6 173.00920 29.384  5.726000 64.40298 1909.250

``` r
# The custom.featureTable() function outputs a dataframe formatted featuretable as well as an MSnbase object.
# The "level" column shows the level of each feature.
featureTable <- custom.featureTable(dir = dir, FTdir = FTdir)
head(featureTable)
```

    ##          mz       rt    rtmin    rtmax     maxo sample level
    ## F1 44.99859 1444.718 1435.434 1469.897 82181.38      3    PP 
    ## F2 44.99863 1272.449 1211.726 1281.954 77924.63      3    PP
    ## F3 44.99866 1680.282 1674.256 1693.137 43457.13      1    PP
    ## F4 44.99866 1680.282 1674.256 1693.137 43457.13      1    PP
    ## F5 44.99868 1274.340 1201.896 1290.801 80583.88      1    PP
    ## F6 44.99868 1274.340 1201.896 1290.801 80583.88      1    PP

# Part 3: MS2 recognition Feature Identification

After PP features have been extracted, MR feature
identification can be performed. This step is optional.

``` r
# Find MR features and add them to the original feature table.
featureTable <- find.level3features(data = MSdata)
```

# Part 4: Add Target Features

After extracting all PP and MR features, users may also want to
extract features from an in-house library of standard metabolites with
given mass and retention time information. If these features already
exists in the original feature table as either PP or MR
features, they will not be added. Otherwise, these features will be
explored from the raw data and added to the feature table with “target”
as their level if exist. The in-house metabolites library must be in CSV
format with ONLY the columns “mz” and “rt” in order. This step is
optional.

``` r
# Directory containing the in-house standard metabolite library.
tarFTdir = "X:/Users/JPAtest_20210330"
# Name of the in-house standard metabolite library in CSV format.
tarFTname = "LibraryCSVHILIC-.csv"

# Sample format of in-house standard metabolite library.
setwd(tarFTdir)
head(read.csv(tarFTname, header = T, stringsAsFactors = F))
```

    ##         mz  rt
    ## 1  74.0237 642
    ## 2  87.0088 708
    ## 3  88.0393 444
    ## 4  88.0393 444
    ## 5  89.0233 468
    ## 6 101.0244 720

``` r
# Extract target feature and add them to the original featuretable if they do no exist already.
featureTable <- find.targetfeatures(data = MSdata, tarFTdir = tarFTdir, tarFTname = tarFTname)
# Target features that are newly added.
head(featureTable[featureTable$level == "target",])
```

    ##          mz  rt rtmin rtmax maxo sample  level
    ## F68 88.0393 444   444   444  402      1 target
    ## F69 88.0393 444   444   444  402      1 target
    ## F70 88.0393 444   444   444  330      2 target
    ## F71 88.0393 444   444   444  330      2 target
    ## F72 88.0393 444   444   444  512      3 target
    ## F73 88.0393 444   444   444  512      3 target

# Part 5: Sample Alignment

For multi-sample analysis, sample alignment is performed after all level
features have been identified. Do not perform sample
alignment for single-sample analysis.

``` r
# Update feature table to aligned feature table for multi-sample analysis.
featureTable <- peak.alignment(data = MSdata, bw = 5, minfrac = 0.5, mzwid = 0.015,
                               minsamp = 1, max = 100, quantitative.method = "maxo")
head(featureTable)
```

    ##          mz         rt    rtmin      rtmax NISTplasmaDDARP(-)1.mzXML
    ## F1 44.99871 1428.82483 1428.825 1440.01100                 82195.130
    ## F2 44.99873 1680.73462 1670.793 1682.39300                 43457.130
    ## F3 56.99601 1688.53400 1678.229 1690.69250                  4393.375
    ## F4 56.99611   59.37352   54.225   64.03593                  5108.125
    ## F5 60.99350 1686.36900 1678.229 1687.27600                  2136.125
    ## F6 61.98857 1239.31238 1238.709 1239.31238                 11225.500
    ##    NISTplasmaDDARP(-)2.mzXML NISTplasmaDDARP(-)3.mzXML
    ## F1                 65772.250                     0.000
    ## F2                 39078.130                 49169.380
    ## F3                  4224.875                  4341.875
    ## F4                  5882.500                  5351.500
    ## F5                  2079.500                  1892.375
    ## F6                     0.000                 10942.250

# Part 6: MS2 Annotation

Feature annotation is performed by two functions: “ms2.tofeaturetable”
for assigning MS2 to each corresponding feature, and then
“feature.annotation” for comparing each MS2 against a standard MS2
library for spectra similarity based metabolite annotation. The standard
library used to perform annotation must be in msp format.

``` r
# The ms2.tofeaturetable() function assigns MS2 spectra to the MS1 feature table. It returns a new feature table with additional columns containing MS2 fragment information.
featureTable <- ms2.tofeaturetable(data = MSdata, featureTable = featureTable)

# Now, use the feature.annotation() function to annotate features in the feature table against a standard database in msp format. This functions returns a feature table containing additional columns with annotation information.
lib_directory <- "X:/Users/Library" # directory containing the library file
lib_name <- "convertedLibraryNeg.msp" # name of the library file
featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.1)
# Display some annotated features.
head(featureTable[featureTable$Annotation != "unknown",])
```

    ##            mz         rt      rtmin      rtmax NISTplasmaDDARP(-)1.mzXML
    ## F22  96.96042   33.25266   32.06100   34.44432                     0.000
    ## F23  96.96973 1464.73007 1462.40814 1467.05200                  1568.000
    ## F38 112.98560 1439.97046 1431.12700 1443.37700                  6372.125
    ## F60 129.05550   75.95423   70.53881   79.41445                 10862.880
    ## F61 130.08775   48.60897   44.05600   50.92079                 25484.000
    ## F70 137.02067   92.47435   91.01100   93.93770                 10328.000
    ##     NISTplasmaDDARP(-)2.mzXML NISTplasmaDDARP(-)3.mzXML MS2_match
    ## F22                 20666.000                  27628.00      TRUE
    ## F23                  1624.125                      0.00      TRUE
    ## F38                  6304.625                      0.00      TRUE
    ## F60                 11561.880                  10012.13      TRUE
    ## F61                 19459.250                  15849.88      TRUE
    ## F70                  9604.000                      0.00      TRUE
    ##                               MS2mz             MS2int PeaksCount fromFile
    ## F22 60.9939;78.9591;79.9571;96.9602 838;1960;5378;6714          4        3
    ## F23                 76.9707;78.9595            816;836          2        2
    ## F38                  44.9989;68.996          1608;1328          2        3
    ## F60                          44.999                900          1        2
    ## F61                 44.999;130.0873           436;1650          2        1
    ## F70  65.041;76.9703;93.0346;94.9805 900;3602;10898;816          4        3
    ##                     Annotation   DPscore
    ## F22            Phosphoric acid 0.9701192
    ## F23            Phosphoric acid 0.1128625
    ## F38       Trifluoroacetic acid 0.4502731
    ## F60 3-Methyl-2-oxovaleric acid 1.0000000
    ## F61  Unknown (carbon number 6) 0.9668160
    ## F70             Salicylic acid 0.9462781

# Part 7: Export EIC

EICs can be exported for all level peak picking, MS2 recognition, target, or aligned features
individually. For multi-sample analysis, if users want to export different levels of features’ EIC separately, then this step must be
performed before sample alignment. Otherwise, only aligned features’ EIC
can be exported because the information of feature levels will
not remain in the aligned feature table.

``` r
# Directory to export EICs to.
plotdir <- "X:/Users/JPAtest_20210330"

# Draw and export EICs of all level PP features with MS2.
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 1, smooth = 2)
# Draw and export EICs of all level PP features without MS2.
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 2, smooth = 2)
# Draw and export EICs of all level MR features.
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = 3, smooth = 2)
# Draw and export EICs of all level target features.
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = "target", smooth = 2)

# ONLY for multi-sample analysis AND after alignment, draw and export EICs of all aligned features.
plot.features(dir = plotdir, featureTable = featureTable, data = MSdata, plot.type = "multi", smooth = 2)
```

# Part 8: CAMERA Annotation

At the end of the JPA workflow, users can use CAMERA to identify adduct
and isotope relationships in the feature table. Note: CAMERA annotation
can ONLY be used for single-sample `JPA` analysis.

``` r
# Perform CAMERA annotation.
featureTable <- adduct.isotope.annotation(featureTable = featureTable, data = MSdata, polarity = "negative")
```

# Part 9: Additional Details and Notes
