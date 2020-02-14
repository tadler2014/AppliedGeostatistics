# Cross-correlations,	cross-covariances	and	cross-semivariance

# Created by: Thomas Adler
# Last Updated: September 17th, 2019

# The following script utilizes the "Berea Sandstone" data set to caluclate cross-statistics for
# a	given	lag	distance.	The	Berea	sandstone	is	a	section	of	rock	where	resistivity	
# and	permeability	were	measured	on	extremely	small	intervals	(~3	mm	increments).

# This script is an extension of 'LagStatistics.R'. 
# Unlike 'LagStatistics.R', this code in this script is designed for cross-statistics.

# Note that this script needs to be run twice for the horizontal and vertical directions by
# changing the di (line 39) and dc (line 43) values 


#-------------------- Data Input------------------------
library(readxl)

# define file location
myFile <- "~/Desktop/R/CE269 Applied Geostatistics/Homework 2/Burea_subset.xls"

# read data file
myData <- read_excel(myFile)
# view dataset structure
str(myData)

#-------------------- Define Input Parameters for Lag Calculations -----------------------
# Define grid system
x <- myData$`x (mm)`
y <- myData$`y (mm)`

# Define primary variable of interest
v <- myData$`Permeability (mD)`
# Define secondary variable of interest
v1 <- myData$`Water Resisitivity ohm-m`
# Define direction of interest (di) 
#   di <- x for horizontal assessment
#   di <- y for vertical assessmemt
di <- x
# Define direction to keep constant (dc)
#   dc <- y for horizontal assessment
#   dc <- x for vertical assessmemt
dc <- y

# Define lags
lag = c(1:90)

#-------------------- Create Data Frame for Desired Lag Distances -----------------------
#Create Blank Data Frame of Combinations 
ldf <- data.frame(matrix(ncol = 6, nrow = 0))
n <- c("dc","dlow", "dhigh","vlow", "v1high", "dist")
colnames(ldf) <- n
#Parameter definitions
# dlow <- coordinate position 1
# dhigh <- coordinate position 2
# vlow <- variable of interest at position 1
# v1high <- secondary variable of interest at position 2
# dist <- distance (or lag) between positions 1 and 2

#Fill Blank Data Frame 
dcu <- unique(dc)
j=1
for (k in 1:length(lag)){
  #Create dataset of 
  for (i in 1:length(dcu)){
    # take sample of data set
    sub <- data.frame("dir interest" = di, "dir constant" = dc, "v" = v, "v1"=v1)
    # subset the data by keeping one dimension constant
    sub <- subset(sub, dc == dcu[i])
    if (length(sub$v) > 1){
      # Every combination of locations and variables 
      df1 <- t(combn(sub$dir.interest,2))
      df2 <- t(combn(sub$v,2))
      df3 <- t(combn(sub$v1,2))
      df <- cbind(df1,df2,df3)
      df <- data.frame("dc"=dcu[i],"d low" = df[,1], "d high" = df[,2], "v low" = df[,3], "v1 high" = df[,6])
      # distance between every combination of locations
      df$dist <- abs(df$d.low-df$d.high)
      # subset the data by distances that equal the desired lag
      sdf <-subset(df, dist == lag[k])
      if (!dim(sdf)[1] == 0) {
        # add data to the lag data frame (ldf)
        l<-dim(sdf)[1]
        ldf[j:(j+l-1),] <- sdf
        j = j+l
      }
    }
  }
}

#-------------------- Statistical Assessment of Lag Distances -----------------------
library(dplyr)

# Make sure that the pearson, covariance and semivariance functions are called in.

# Split data by lag
sp <- split(ldf, ldf$dist)

# Calculate Pearson Correlation (p)
p <- lapply(names(sp), function(x) pearson(sp[[x]][["vlow"]],sp[[x]][["v1high"]]))
p <- unlist(p)

# Calculate Covariance (cov)
cov <- lapply(names(sp), function(x) covariance(sp[[x]][["vlow"]],sp[[x]][["v1high"]]))
cov <- unlist(cov)

# Calculate semivariance (s)
s <- lapply(names(sp), function(x) semivariance(sp[[x]][["vlow"]],sp[[x]][["v1high"]]))
s <- unlist(s)

# Combining result
results <- data.frame("Lag"= lag,"p"= p,"cov"= cov,"s"= s)
results <- round(results, digits=3)

# Export Results as an Excel Table
# write.csv(results, file=paste("HW3Results_Horizontal.csv"))

#-------------------- Plotting Results -----------------------
# Plotting statistical parameters against lag distance
par(mfrow=c(1,3))
plot(results$Lag,results$p,ylab="Pearson's Correlation Coeff (r)",xlab="Lag (mm)",col = "black",pch = 16,cex=1,type="o")
plot(results$Lag,results$cov,ylab="Covariance (cov)",xlab="Lag (mm)",col = "black",pch = 16,cex=1,type="o")
plot(results$Lag,results$s,ylab="Semivariance",xlab="Lag (mm)",col = "black",pch = 16,cex=1,type="o")

# Plotting h-scatterplots at different lag steps
names(sp) <- c("3 mm", "6 mm","9 mm","12 mm","15 mm","18 mm")
par(mfrow=c(3,2))
lapply(1:length(lag), function(x) {plot(sp[[x]][["vlow"]],sp[[x]][["v1high"]],
                                       xlab="Permeability (x)",ylab="Resistivity (x+h)",
                                       xlim=c(190,592),ylim=c(240,260),main=names(sp[x]),
                                       abline(lm(sp[[x]][["v1high"]] ~ sp[[x]][["vlow"]])))
       text(570,259,label=bquote(r==.(results$p[x])))
       text(570,257,label=bquote(cov==.(results$cov[x])))
       text(570,255,label=bquote(s==.(results$s[x])))})
       
                                    
