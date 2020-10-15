
rm(list=ls()) #clear the workspace
setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes")
data<-read.csv("vent_fields_all.csv")
library("spatstat")
library("proj4")
source("create_landscape_func.R")

unique(data$Activity)
#if we wanted to subset for only active vents
#data<-data[data$Activity == "active, confirmed" | data$Activity == "active, inferred",]

#remove any rows with missing coordinates
data[!is.na(data$Latitude) | !is.na(data$Longitude),]

unique(data$Region)
EPR.data<-data[(data$Region == "N EPR" | data$Region == "S EPR"),]
MAR.data<-data[(data$Region == "N MAR" | data$Region == "S MAR"),]

EPR.coordinates<-data.frame(EPR.data$Longitude, EPR.data$Latitude)
EPR.interpatch.distances<-dist(EPR.coordinates, method="euclidean")
hist(EPR.interpatch.distances)
plot(EPR.data$Longitude, EPR.data$Latitude)

################
install.packages("rgdal")
library(ggplot2)
#install.packages("ggmap")
#library("ggmap")
#library("tidyverse")
#map.world<-get_map("Tokyo")
#?register_google

#data(state)
#s <- project(state.center, "+proj=merc", inverse=TRUE)
#project(EPR.coordinates, proj, inverse = TRUE, degrees = TRUE, silent = FALSE,
#        ellps.default="sphere")
########################


#print(ggplot(EPR.coordinates, aes(Longitude,Latitude)) + theme_classic() + geom_point()) #+ geom_point(aes(size = areas)))
#calculate distribution of nndists for that landscape
EPR.nndists<-nndist(EPR.data$Longitude, EPR.data$Latitude)
hist(EPR.nndists)
summary(EPR.nndists)

MAR.coordinates<-data.frame(MAR.data$Longitude, MAR.data$Latitude)
MAR.interpatch.distances<-dist(MAR.coordinates, method="euclidean")
hist(MAR.interpatch.distances)
#calculate distribution of nndists for that landscape
MAR.nndists<-nndist(MAR.data$Longitude, MAR.data$Latitude)
hist(MAR.nndists)
summary(MAR.nndists)


landscape<-create.landscape(n.patches=50, landscape.type="linear", landscape.limit=20, 
                            patch.distribution="clustered", areas.distribution="lognormal", 
                            areas.limit=2, clustering.iters=100)
head(landscape)
print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
landscape.coordinates<-data.frame(landscape$x.coord,landscape$y.coord)
landscape.interpatch.distances<-dist(landscape.coordinates, method="euclidean")
hist(landscape.interpatch.distances)
#calculate distribution of nndists for that landscape
landscape.nndists<-nndist(landscape$x.coord, landscape$y.coord)
hist(landscape.nndists)
summary(landscape.nndists)

