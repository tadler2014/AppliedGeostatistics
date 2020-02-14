# Contour Map of Variogram Values

# See page 151 of Isaaks et. al for an example

# Created by: Thomas Adler
# Last Updated: October 29th, 2019

#-------------Data Input -----------------
myData <- read.csv("berea_subset.csv")

#-------------Define Semivariogram Parameters-----------------
#Define variables of interest
x <- myData$x_mm # Mesh direction 1
y <- myData$y_mm# Mesh direction 2
v <- myData$Permeability_mD # Variable of interest

#-------------Equal Split of Data Semivariograms -----------------
# Create a semivariogram with an	equal	number	of	data	points	in	each	bin.
#Example: Semivariogram w/ 100 bins with equal sample sizes
myData2 <- semivarES(x,y,v,100)

#-------------Contour Map Bins -----------------
# Pull data
myData3 <- myData2[["AllData"]]
myData3 <- subset(myData3, select=c("s1", "xdist1","ydist1"))
rm(myData2)


# define interval breaks by selected bin width
gridsquares_x <- seq(from = -300, to = 0, by = 20)
gridsquares_ypos <- seq(from = 0, to = 90, by = 10)
gridsquares_yneg <- seq(from = -90, to = 0, by = 10)

gridsquares_y <- seq(from = -90, to = 90, by = 20)


#-------------Contour Map -----------------
# Calculate two corners of the contour map
contourData <- data.frame("svar"=numeric(length(gridsquares_x)*length(gridsquares_y)),"x"=numeric(length(gridsquares_x)*length(gridsquares_y)),"y"=numeric(length(gridsquares_x)*length(gridsquares_y)))
k=1
for (i in 1:length(gridsquares_x)){
  temp <- subset(myData3, gridsquares_x[i] < myData3$xdist1 & myData3$xdist1 < gridsquares_x[i+1])
  for (j in 1:length(gridsquares_y)){
    temp1 <- subset(temp, gridsquares_y[j] < temp$ydist1 & temp$ydist1 < gridsquares_y[j+1])
    contourData$x[k] <- (gridsquares_x[i]+gridsquares_x[i+1])/2
    contourData$y[k] <- (gridsquares_y[j]+gridsquares_y[j+1])/2
    contourData$svar[k] <- mean(temp1$s1)
    k=k+1
  }
}
rm(temp,temp1)
contourData <- na.omit(contourData) 
# Calculate the other two corners of the contour
contourData2 <- contourData
contourData2$x <- contourData2$x*-1
contourData2$y <- contourData2$y*-1
# Bind data
contourData3 <- rbind(contourData, contourData2)
# Plot Results
require(lattice)
contourplot(svar ~ x * y, data = contourData3,cuts = 20,xlim=c(-300,300),ylim=c(-90,90),asp=0.25)




