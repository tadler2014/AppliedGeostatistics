# Equal Split semivariogram calculator

# Function: Creates a semivariogram with an	equal	number	of	data	points	in	each	bin.

# Created by: Thomas Adler
# Last edited: September 29th, 2019

semivarES <- function(x,y,v,bins){
  # Step 1: Calculate euclidean distance (d) between each point in the x-y matrix
  d = array(0, dim=c(length(x),length(y))) #Distance Matrices (square)
  for (i in 1:length(x)){  #Location 1
    for (j in 1:length(y)){ #Location 2
      if (j==i) { 
        # Ignore distances of 0
        d[i,j] <- NA 
      } else {
        # Euclidean Distance
        d[i,j] <- sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2 ) 
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
        s[i,j] <- 0.5*(v[j]-v[i])^2 
      }}}
  s[upper.tri(s)] <- NA #Remove replicates
  #Vectorize Matrix
  s1 <- as.vector(s)
  s1<- s1[!is.na(s1)] #vectorize matrix
  
  # Step 3: Create data frame with euclidean distance (d1) and semivariance (s1)
  myData2 <- data.frame(d1, s1)
  myData2 <- myData2[order(d1),] #order the dataframe
  
  
  # Step 4: Find mean at equal intervals equal interval bins
  # split data into equal interval bins
  d2 <- split(myData2$d1, ceiling(seq_along(myData2$d1)/(length(myData2$d1)/bins)))
  s2 <- split(myData2$s1, ceiling(seq_along(myData2$d1)/(length(myData2$d1)/bins)))
  # find mean
  d3 <- sapply(d2, mean)
  s3 <- sapply(s2, mean)
  
  # Step 5: Plot results and store data
  #plot(d1,s1,xlab="Distance(mm)",ylab="Variance(mD)")
  #xvals=vals[["x"]]
  #abline(v=xvals, col=c("blue"), lty=1, lwd=c(1))
  plot(d1,s1,xlab="Distance(mm)",ylab="Variance(mD)")
  plot(d3,s3,xlab="Distance(mm)",ylab="Variance(mD)")
  
  # store data
  mysemiData <- data.frame("Distance"=d3, "Semivariance"=s3)
  return(mysemiData)
  rm (d, d1, d2, d3, myData2, s, s1, s2, s3, i, j, N)
}
