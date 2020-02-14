# Spatial Statistics Functions

# Created by: Thomas Adler
# Last Edited: October 11th, 2019

# Correlation
pearson <- function(x,y){
  r <- sum((x-mean(x))*(y-mean(y)))/sqrt(sum((x-mean(x))^2)*sum((y-mean(y))^2))
}

# Covariance
covariance <- function(x,y){
  cov <- sum((x-mean(x))*(y-mean(y)))/(length(x)-1)
}

# Semivariance
semivariance <- function(x,y){
  s <- (1/(2*length(x)))*sum((x-y)^2)
}

# Cross-Correlation
crosspearson <- function(x,y){
  r <- sum((x-mean(x))*(y-mean(y)))/sqrt(sum((x-mean(x))^2)*sum((y-mean(y))^2))
}

# Cross-Covariance
crosscovariance <- function(x,y){
  cov <- sum((x-mean(x))*(y-mean(y)))/(length(x)-1)
}

# Cross-Semivariance
crosssemivariance <- function(x,y){
  s <- (1/(2*length(x)))*sum((x-y)^2)
}




