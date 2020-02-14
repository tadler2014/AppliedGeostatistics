# Semivariogram Equations

#--------------------Equal Split semivariogram calculator----------------
# Function: Creates a semivariogram with an	equal	number	of	data	points	in	each	bin.

# Created by: Thomas Adler
# Last edited: September 29th, 2019

semivarES <- function(x,y,v,bins,angle){
  # Step 1: Calculate euclidean distance (d) between each point in the x-y matrix
  d = array(0, dim=c(length(x),length(y))) #Distance Matrices (square)
  theta = array(0, dim=c(length(x),length(y))) #Angle Matrices (square)
  xdist = array(0, dim=c(length(x),length(y))) #Angle Matrices (square)
  ydist = array(0, dim=c(length(x),length(y))) #Angle Matrices (square)
  for (i in 1:length(x)){  #Location 1
    for (j in 1:length(y)){ #Location 2
      if (j==i) { 
        # Ignore distances of 0
        d[i,j] <- NA 
        theta[i,j] <- NA 
        xdist[i,j] <- NA 
        ydist[i,j] <- NA 
      } else {
        # Euclidean Distance
        d[i,j] <- sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2 )
        # angle between points 
        theta[i,j] = (acos((x[j]-x[i])/sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2))*180)/(pi)
        # x distance
        xdist[i,j] <- (x[j]-x[i])
        # y distance
        ydist[i,j] <- (y[j]-y[i])
      }}}
  d[upper.tri(d)] <- NA #Remove replicates
  theta[upper.tri(theta)] <- NA #Remove replicates
  xdist[upper.tri(xdist)] <- NA #Remove replicates
  ydist[upper.tri(ydist)] <- NA #Remove replicates
  #Vectorize matrix
  d1 <- as.vector(d)
  d1<- d1[!is.na(d1)] 
  theta1 <- as.vector(theta)
  theta1<- theta1[!is.na(theta1)] 
  xdist1 <- as.vector(xdist)
  xdist1<- xdist1[!is.na(xdist1)] 
  ydist1 <- as.vector(ydist)
  ydist1<- ydist1[!is.na(ydist1)] 
  
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
  
  # Step 3: Create data frame with euclidean distance (d1), angle (theta1) and semivariance (s1)
  myData2 <- data.frame(d1, theta1, s1, xdist1, ydist1)
  myData2 <- myData2[order(d1),] #order the dataframe
  
  # Step 3a: if angle is selected, take a subset
  if (!missing(angle)){
    myData2 <- subset(myData2, theta1 == angle)
  }
  
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
  #plot(d1,s1,xlab="Distance(mm)",ylab="Variance(mD)")
  plot(d3,s3,xlab="Distance(mm)",ylab="Variance(mD)")
  
  # store data
  binnedsemiData <- data.frame("Distance"=d3, "Semivariance"=s3)
  allsemiData <- myData2
  
  mysemiData <- list("BinnedData" = binnedsemiData,"AllData" = allsemiData)
  
  return(mysemiData)
  rm (d, d1, d2, d3, myData2, s, s1, s2, s3, i, j, N)
}

#--------------------User Defined semivariogram calculator---------------
# Function: Allows user to select bins for semivariograms

semivarUB <- function(x,y,v,n_user,angle){
  # Step 1: Calculate euclidean distance (d) between each point in the x-y matrix
  d = array(0, dim=c(length(x),length(y))) #Distance Matrices (square)
  theta = array(0, dim=c(length(x),length(y))) #Angle Matrices (square)
  for (i in 1:length(x)){  #Location 1
    for (j in 1:length(y)){ #Location 2
      if (j==i) { 
        # Ignore distances of 0
        d[i,j] <- NA 
        theta[i,j] <- NA 
      } else {
        # Euclidean Distance
        d[i,j] <- sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2 )
        # angle between points 
        theta[i,j] = (acos((x[j]-x[i])/sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2))*180)/(pi)
      }}}
  d[upper.tri(d)] <- NA #Remove replicates
  theta[upper.tri(theta)] <- NA #Remove replicates
  #Vectorize matrix
  d1 <- as.vector(d)
  d1<- d1[!is.na(d1)] 
  theta1 <- as.vector(theta)
  theta1<- theta1[!is.na(theta1)] 
  
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
  
  # Step 3a: if angle is selected, take a subset
  if (!missing(angle)){
    myData2 <- subset(myData2, theta1 == angle)
  }
  
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


#--------------------Bin Width semivariogram calculator------------------
# Function: Creates a semivariogram with user defined bin widths

semivarEW <- function(x,y,v,bin_width,angle){
  # Step 1: Calculate euclidean distance (d) between each point in the x-y matrix
  d = array(0, dim=c(length(x),length(y))) #Distance Matrices (square)
  theta = array(0, dim=c(length(x),length(y))) #Angle Matrices (square)
  for (i in 1:length(x)){  #Location 1
    for (j in 1:length(y)){ #Location 2
      if (j==i) { 
        # Ignore distances of 0
        d[i,j] <- NA 
        theta[i,j] <- NA 
      } else {
        # Euclidean Distance
        d[i,j] <- sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2 )
        # angle between points 
        theta[i,j] = (acos((x[j]-x[i])/sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2))*180)/(pi)
      }}}
  d[upper.tri(d)] <- NA #Remove replicates
  theta[upper.tri(theta)] <- NA #Remove replicates
  #Vectorize matrix
  d1 <- as.vector(d)
  d1<- d1[!is.na(d1)] 
  theta1 <- as.vector(theta)
  theta1<- theta1[!is.na(theta1)] 
  
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
  
  # Step 3a: if angle is selected, take a subset
  if (!missing(angle)){
    myData2 <- subset(myData2, theta1 == angle)
  }
  
  # Step 4: Find mean at defined bin widths
  # define interval breaks by selected bin width
  b <- seq(from = 0, to = max(x), by = bin_width)
  # add max value of lag distance to breaks if the max value isn't already included
  if (tail(b, n=1) != max(x)){
    b <- union(b, max(x))
  }
  # split data into bins by bin width
  d2 <- split(myData2$d1, cut(myData2$d1, breaks = b))
  s2 <- split(myData2$s1, cut(myData2$d1, breaks = b))
  # find mean
  d3 <- sapply(d2, mean)
  s3 <- sapply(s2, mean)
  
  # Step 5: Plot Results and store data
  # plot results
  #plot(d1,s1,xlab="Distance(mm)",ylab="Variance(mD)")
  plot(d3,s3,xlab="Distance(mm)",ylab="Variance(mD)")
  
  # store data
  mysemiData <- data.frame("Distance"=d3, "Semivariance"=s3)
  return(mysemiData)
  rm (d, d1, d2, d3, myData2, s, s1, s2, s3, i, j, N, b)
}

