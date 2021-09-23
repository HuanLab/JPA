#noise average and sd estimation#
directory <- "F:/Jian_Guo"
setwd(directory)
filename <- list.files(pattern = ".mzXML")
ms1data <- readMSData(files = filename[q], mode = "onDisk", msLevel. = 1)
mzData <- mz(ms1data)
intData <- intensity(ms1data)
currmzRange <- c(mzmin, mzmax)
    tmpMZdata <- mzData
    tmpINTdata <- intData
    for(j in 1:length(mzData)){
      index <- which(tmpMZdata[[j]] >= currmzRange[1] & tmpMZdata[[j]] < currmzRange[2])
      tmpMZdata[[j]] <- tmpMZdata[[j]][index]
      tmpINTdata[[j]] <- tmpINTdata[[j]][index]
    }
    
    eicINTraw <- c()
    eicINT <- c()
    eicRT <- c()
    for(k in 1:length(mzData)){
      if(length(tmpINTdata[[k]]) > 0){
        eicINTraw[k] <- mean(tmpINTdata[[k]])
      }else{
        eicINTraw[k] <- 0
      }
      eicRT[k] <- rtime[k]
    }
    if(sum(eicINTraw != 0) == 0) next()
    eicINT <- peak_smooth(eicINTraw)
    eicNon0 <- sort(eicINT[eicINT > 0])
    if(length(eicNon0) > 10){
      for(x in seq(10,length(eicNon0), 10)){
        sd <- sd(eicNon0[1:x])
        blk <- sum(eicNon0[1:x])/x
        thres <- blk + 3*sd
        if(x+1 <= length(eicNon0)){
          if(eicNon0[x+1] >= thres) break()
        }
      }
    }
signal <- blk
variance <- sd


#the noise threshold determination model#
hz <- 5
signal <- 100
variance <- 40
for (x in 1:10) {
 met.threshold <- rep(FALSE, 100)
 for(i in 1:length(met.threshold)) {
  d <- rnorm(hz*60, mean = signal, sd = variance)
  d
  if(max(d)>x*signal) {met.threshold[i] <- TRUE}
 }
 if (length(met.threshold[met.threshold == FALSE]) > 95){
 break
 }
}
print(paste("use threshold", x))
