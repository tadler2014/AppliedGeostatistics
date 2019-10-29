# Ordinary Kriging Function

# Created by: Thomas Adler
# Last Updated: October 8th, 2019

# Inputs:
# myData <- values at known locations(x,y,v1)
# grid <- expanded grid of interest

OKriging <- function(myData,grid,myModel){
  # Step 0: Create empty matrix for final results
  filled.grid <- data.frame("x"=numeric(),"y"=numeric(),"var1.predict"=numeric(),"var1.var"=numeric())
  # Step 1: Semivariance Matrix (C)
  # Calculate euclidean distance (d) between each known point in the x-y matrix
  d = array(0, dim=c(length(myData$x),length(myData$y))) #Distance Matrices (square)
  for (i in 1:length(myData$x)){  #Location 1
    for (j in 1:length(myData$y)){ #Location 2
      # Euclidean Distance
      d[i,j] <- sqrt((myData$x[j]-myData$x[i])^2 + (myData$y[j]-myData$y[i])^2 ) 
    }}
  # Calculate semivariance between each known point in the x-y matrix (the C matrix)
  C = array(0, dim=c(nrow(d)+1,ncol(d)+1)) #Distance Matrices (square)
  for (i in 1:length(myData$x)){  #Location 1
    for (j in 1:length(myData$y)){ #Location 2
      # Semivariance
      C[i,j] <- myModel(d[i,j]) 
    }}
  # Finalize C matrix
  C[nrow(d)+1,]=1
  C[,ncol(d)+1]=1
  C[nrow(d)+1,ncol(d)+1]=0
  
  # Step 2: Distance Vector (D)
  for (k in 1:nrow(grid)){
    # Define unknown location
    x0 = grid$x[k]
    y0 = grid$y[k]
    # If at location of known points, fill in and move along
    if (any(myData[, "x"] == x0 & myData[, "y"] == y0)){
      point <- which((myData$x == x0) & (myData$y == y0))
      filled.grid[k,1] <- x0
      filled.grid[k,2] <- y0
      filled.grid[k,3] <- myData$v1[point]
      filled.grid[k,4] <- 0
    } else {
      # Calculate euclidean distance (d) between each known point and the unknown.
      d1 = vector(mode="numeric")
      for (i in 1:length(myData$x)){  #Location 1
        d1[i] <- sqrt((x0-myData$x[i])^2 + (y0-myData$y[i])^2) 
      }
      # Calculate semivariance between each known point and the unknown.
      D = vector(mode="numeric")
      for (i in 1:length(d1)){  #Location 1
        D[i] <- myModel(d1[i])
      }
      # Finalize D matrix
      D[length(d1)+1]=1
      
      # Step 3: Weight Vector (W)
      # Take the inverse of C
      invC <- solve(C)
      # Multiply inverse of C by D
      w = invC %*% D
      w[length(w)] = abs(w[length(w)])
      
      # Step 4: Estimate of value (v0) and error variance (v0_var) at select point
      filled.grid[k,1] <- x0
      filled.grid[k,2] <- y0
      filled.grid[k,3] = sum(w[1:(length(w))-1]*myData$v)
      filled.grid[k,4] = sum(w[1:(length(w))-1]*D[1:(length(D))-1])+w[length(w)]
    }
  }
  # Step 5: Plot Results
  # Predicted Variables
  matrix.predict <- xtabs(var1.predict~x+y, data=filled.grid)
  filled.contour(x,y,matrix.predict,color.palette = terrain.colors,plot.axes = { 
    points(x = myData$x, y = myData$y)},
    key.title = title(main="   Na+ [mg/L]",cex.main=0.70))
  
  # Variance
  matrix.var <- xtabs(var1.var~x+y, data=filled.grid)
  filled.contour(x,y,matrix.var,color.palette = heat.colors,plot.axes = { 
    points(x = myData$x, y = myData$y)},
    key.title = title(main="   Error Variance",cex.main=0.70))
  return(filled.grid)
}
