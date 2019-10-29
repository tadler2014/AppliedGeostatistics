# Homework 5 Driver

# Created by: Thomas Adler
# Last Updated: October 8th, 2019

# Code should operate for both semivariance and covariance
#-------------Data Input -----------------
# Part 1: Create Dataframe of Known Data
myData <- data.frame("x"=c(5,10,18,3,7),"y"=c(5,7,13,15,2),"v1"=c(150,225,400,100,800))
# Note: For the ordinary kriging function to work it is important that the known data be
#       organized in a data frame like such:
# myData$x <- x-direction
# myData$y <- y-direction
# myData$v1 <- variable of interest

# Part 2: Create Grid
# Define mesh 
meshsize = 1
# Define 2-D bounds and expand grid
x = seq(0,20,by=meshsize)
y = seq(0,20,by=meshsize)
grid <- expand.grid(x,y)
names(grid) <- c("x", "y")

# Part 3: Create Semivariogram Model
# Model Options
exponential = function(h){c0+c*(1-exp(-h/a))} #Equation
# Model parameters
c0 = 50 #Nugget 0 
w = 15 #Sill 15
a = 50 #Range 10
c = w-c0

#-------------Ordinary Kriging-----------------
filledgrid <- OKriging(myData,grid,exponential)


#-------------Plotting Points with unknown at (7,14)-----------------

plot(y~x,myData,ylim=c(0,20),xlim=c(0,20))
text(y~x,labels=1:5,data=myData, cex=0.9, font=2,pos=4)
par(new=TRUE)
plot(7,14,ylim=c(0,20),xlim=c(0,20),ylab = "y",xlab = "x")
text(7,14,labels=paste("unknown (7,14)"),cex=0.9, font=2,pos=4)

