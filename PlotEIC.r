library(xcms)
directory <- "F:/Jian_Guo/IPA_paper_ver2_20210421/Response_20210805/EICexampleformetaboliteSTD"
setwd(directory)
input.file <- list.files(pattern = ".mzXML")
csv.file <- list.files(pattern = ".csv")
featureTable <- read.csv(csv.file[1], header = T, stringsAsFactors = F)
xraw <- xcmsRaw(input.file[1],profstep=0, mslevel = 1)
plot.matrix <- featureTable
if(nrow(plot.matrix)!=0){
  for(k in 1:nrow(plot.matrix)){
    rt.lower.limit <- plot.matrix$rt[k] - 60
    rt.upper.limit <- plot.matrix$rt[k] + 60
    mass.lower.limit <- plot.matrix$mz[k] - 0.01
    mass.upper.limit <- plot.matrix$mz[k] + 0.01
    mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
    RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
    eeic <- rawEIC(xraw, mzrange=mzRange, rtrange=RTRange) #extracted EIC object
    points <- cbind(xraw@scantime[eeic$scan], eeic$intensity)
    png(file = paste0(rownames(plot.matrix)[k], ".png"),
        width = 480, height = 480)
    eic <- plot(points, type="l", main= plot.matrix$name[k], xlab="Seconds",
                ylab="Intensity", xlim=RTRange)
    dev.off()
  }
}
