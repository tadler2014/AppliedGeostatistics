# User Defined semivariogram calculator

# Created by: Thomas Adler
# Last edited: September 24th, 2019

semivarUB <- function(x,y,v,n_user){
  # Step 1: Calculate euclidean distance (d) between each point in the x-y matrix
  d = array(0, dim=c(length(x),length(y))) #Distance Matrices (square)
  for (i in 1:length(x)){  
    for (j in 1:length(y)){ 
      if (j==i) {
        # Ignore distances of 0
        d[i,j] <- NA
      } else {
        # Euclidean Distance
        d[i,j] <- sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2 ) #Euclidean Distance
      }}}
  d[upper.tri(d)] <- NA #Remove replicates
  #Vectorize matrix
  d1 <- as.vector(d)
  d1<- d1[!is.na(d1)] 
  
  # Step 2: Calculate semivariance (s) between each point in the x-y matrix
  N = length(v)
  s = array(0, dim=c(N,N))
  for (i in 1:N){  
    for (j in 1:N){ 
      if (j==i) {
        # Ignore variables at same point
        s[i,j] <- NA
      } else {
        # Semivariance
        s[i,j] <- 0.5*(v[j]-v[i])^2 #Semivariance
      }}}
  s[upper.tri(s)] <- NA #Remove replicates
  #Vectorize Matrix
  s1 <- as.vector(s)
  s1<- s1[!is.na(s1)] #vectorize matrix
  
  # Step 3: Create data frame with euclidean distance (d1) and semivariance (s1)
  myData2 <- data.frame(d1, s1)
  myData2 <- myData2[order(d1),] #order the dataframe
  
  # Step 4: Define Bins
  #Choose bin locations(NEW)
  print("Plotting full semivariogram...")
  plot(d1,s1, xlab="Distance(mm)",ylab="Variance(mD)")
  Sys.sleep(1)
  print("Choose first bin window by clicking on the semivariogram")
  Sys.sleep(1)
  print("Wait until the blue line appears before choosing your next bin window")
  Sys.sleep(2)
  vals <- list()
  xvals <- list()
  for (i in 1:n_user){
    #Click on the graph where the bins should be placed.
    vals[[i]]=locator(1) 
    #Show bin dividers
    xvals[[i]] <- vals[[i]][["x"]]
    abline(v=xvals[[i]], col=c("blue"), lty=1, lwd=c(1))
  }
  
  # Step 5: Calculate Bin Averages
  print("Calculating bin averages...")
  Sys.sleep(2)
  d2u<-split(myData2$d1, cut(myData2$d1, c(0,xvals,max(d1)), include.lowest=TRUE))
  s2u<-split(myData2$s1, cut(myData2$d1, c(0,xvals,max(d1)), include.lowest=TRUE))
  d3u <- sapply(d2u, mean)
  s3u <- sapply(s2u, mean)
  plot(d3u,s3u,xlab="Distance(mm)",ylab="Variance(mD)")
  #Plot Bin Averages and Return Data Frame
  mysemiData <- data.frame("Distance"=d3u, "Semivariance"=s3u)
  return(mysemiData)
  rm(d, d1, d2u, d3u, s, s1, s2u, s3u, i, j, N, xvals, vals)
}
