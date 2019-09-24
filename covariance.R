# Covariance Function

# Created by: Thomas Adler
# Last Updated: September 11th, 2019


covariance <- function(x,y){
  cov <- sum((x-mean(x))*(y-mean(y)))/(length(x)-1)
}
