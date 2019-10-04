# Homework 4 Driver

# Created by: Thomas Adler
# Last Updated: September 24th, 2019

#-------------Data Input -----------------
myData <- read.csv("berea_subset.csv")

#-------------Define Semivariogram Parameters-----------------
#Define variables of interest
x <- myData$`x (mm)` # Mesh direction 1
y <- myData$`y (mm)`# Mesh direction 2
v <- myData$`Permeability (mD)` # Variable of interest

#-------------Equal Split of Data Semivariograms -----------------
# Create a semivariogram with an	equal	number	of	data	points	in	each	bin.
#Example: Semivariogram w/ 100 bins with equal sample sizes
myData2 <- semivarES(x,y,v,100)

#-------------User-specified bin width Semivariograms -----------------
# Create a semivariogram with bins of equal length
#Example: Semivariogram w/ 5mm bins
myData3 <- semivarEW(x,y,v,5)

#-------------User-Defined Bins Semivariograms -----------------
#Define variables of interest
#Example: Semivariogram w/ 8 user defined bins
myData3 <- semivarUB(x,y,v,8)

#-------------Define model functions -----------------
#Power Equation
power = function(h,g){w*h^g} 
#Spherical
spher = function(h){
  p1 <- c0+c*(1.5*(h[1:a]/a)-0.5*(h[1:a]/a)^3) #d<a
  p2 <- rep(w, times = length(h[(a+1):tail(h,1)])+1) #d>a
  output <- c(p1,p2) #combine
  return(output)
  } 
#Exponential
expon = function(h){c0+c*(1-exp(-h/a))} #Equation
#Gaussian
gauss = function(h){c0+c*(1-exp(-(h/a)^2))} #Equation
#Cubic
cube = function(h){
  p1 <- c0+c*(7*(h[1:a]/a)^2-8.75*(h[1:a]/a)^3+3.5*(h[1:a]/a)^5-0.75*(h[1:a]/a)^7) #d<a
  p2 <- rep(w, times = length(h[(a+1):tail(h,1)])+1) #d>a
  output <- c(p1,p2) #combine
  return(output)  
  }

#-------------Overlay Models on semivariogram -----------------
#Locate Sill and Range
parms=locator(1)
#x = range
#y = Sill

#Define Model Coefficients
c0 = 0 #Nugget
w = parms[["y"]] #Sill
a = parms[["x"]] #Range
c = w-c0
h = 1:max(x) #Lag Distance

#Plot Semivariogram w/ Models
plot(myData3$Distance, myData3$Semivariance,xlab="Distance(mm)",ylab="Variance(mD)", ylim = c(0, 20000))

# power model
lines(h,power(h,0.01),lty=1,col="red",lwd=3)
# spherical model
lines(h,spher(h),lty=2,col="red",lwd=3)
# exponential model
lines(h,expon(h),lty=3,col="green",lwd=3)
# gaussian line
lines(h,gauss(h),lty=4,col="black",lwd=3)
# cubic line
lines(h,cube(h),lty=4,col="orange",lwd=3)

#Add legend
legend("topright",legend=c("Power (lambda = 0.01)", "Spherical","Exponential","Gaussian","Cubic"),
       col=c("red", "blue","green","black","orange"), lty=1:4, cex=0.8,
       box.lty=0)





