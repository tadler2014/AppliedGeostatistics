# Pearson Correlation Coefficient Function

# Created by: Thomas Adler
# Last Updated: September 11th, 2019

pearson <- function(x,y){
  r <- sum((x-mean(x))*(y-mean(y)))/sqrt(sum((x-mean(x))^2)*sum((y-mean(y))^2))
}
