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
