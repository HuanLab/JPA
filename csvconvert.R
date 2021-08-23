setwd("F:/Jian_Guo/IPA_paper_ver2_20210421/test")
filename <- list.files(pattern = ".csv")
for (y in 1:length(filename)) {
  DF <- read.csv(filename[y])
  newDF <- as.data.frame(cbind(DF[,7], DF[,5]*60, DF[,4]*60, DF[,6]*60, DF[,8]))
  colnames(newDF) <- c("mz", "rt", "rtmin", "rtmax", "Height")
  write.csv(newDF, file = paste(filename[y], "converted.csv",  sep = "_"), row.names = FALSE )
}
