# Semivariance Function

# Created by: Thomas Adler
# Last Updated: September 11th, 2019

semivariance <- function(x,y){
  s <- (1/(2*length(x)))*sum((x-y)^2)
}
