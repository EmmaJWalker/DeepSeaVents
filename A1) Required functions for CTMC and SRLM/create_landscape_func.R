#CREATE LANDSCAPE FUNCTION
#################################################################################################################
#creates a landscape of a specified size with a specified number of patches, with a patch distribution
#of specified level of clustering or uniformity, as either a linear array or 2D landscape of patches
#returns a data table containing the patch.ID (patch specific ID number), areas, coordinates, 
#distance.matrix (giving interpatch distances)
#######################################################################################################

#ARGUMENTS
######################################################################################################
#n.patches = number of patches in desired landscape

#landscape.type = specify the type of patch arrangement (accepts either "linear" or "2D")

#landscape.limit = size of disired landscape (max possible x and y coordinate of a patch)

#patch.distribution = specify whether the desired landscape's patches should be uniformly 
#(arranged on a grid in 2D, -note area limit should then be chosen appropiately to ensure 
#no overlap-), randomly, more clustered or more regularly distributed 
#(accepts either uniform, "clustered", "regular" or "random")

#areas.distribution = how patch areas should be distributed in the landscape (accepts "lognormal"
# -meant to best mimic most natural landscapes but avoiding having a large "mainland" patch that
#completely dominates the landscape within a 100x100 unit landscape- ,"uniform" or "random") 

#areas.limit = for "uniform" area distributions sets the area of all patches, 
#for "random" area distributions sets the maximum area of all patches relative to 1

#clustering.iters = number of times clustering/declustering algorithm should be iterated
#higher the value, the more patches are clustered or regularly spaced
#if it is a "clustered" landscape versus a "regular" landscape
#####################################################################################################

#################################################################################
#GENERATE RANDOM LANDSCAPE OF PATCHES
create.landscape<-function(n.patches, landscape.type, landscape.limit, patch.distribution, 
                           areas.distribution, areas.limit, clustering.iters){
  n.patches<-n.patches
  patch.ID<-c(1:n.patches)
  #assign ID numbers to each patch (This is so that we can track which patches are
  #removed throughout destruction of patches from the landscape)
  if (areas.distribution == "lognormal") { #if areas are to be distributed lognormally
    areas<-rlnorm(n.patches, meanlog = 2, sdlog=1 )
    #patch areas lognormally distributed with a mean of 2 and standard deviation of 1
    #this was chosen to best mimic most natural landscapes with lognormally distributed areas
    #but also set so that we don't get one patch taking up almost all the landscape 
    #(acting like a mainland, that every other patch's dynamics depends very heavily on)
    #within a 100x100 unit landscape
  }
  if (areas.distribution == "uniform"){ #if areas are to be distributed uniformly
    areas<-rep( areas.limit, n.patches) #patch areas will all be at areas.limit
  }
  if (areas.distribution == "random"){ #if areas are to be distributed randomly
    areas<-runif(n.patches, min=1, max=areas.limit) #patch areas will be randomly selected from a 
    #uniform distribution between one and the area limit
  }
  

  radii<-sqrt(areas/pi) #calculate radii of hypothetically circular patches
  
  if (patch.distribution!="uniform"){#if not arranging aranging pactches uniformly (grid in 2D)
  
  if (landscape.type == "linear"){ #if creating a linear landscape
    x.coord<-rep(0, n.patches) #patches will be arranged vertically along the y axis (x = 0)
  } else if (landscape.type == "2D"){
    x.coord<-runif(n.patches, min=0, max=landscape.limit) 
    #pick a random number for the x coordinate of each patch between 0 and the extent
    #of the landscape
  }
  y.coord<-runif(n.patches, min=0, max=landscape.limit) 
  #pick a random number for the y coordinate of each patch between 0 and the extent 
  #of the landscape
  coordinates<-data.frame(x.coord,y.coord) 
  #the landscape of patch locations is given by the chosen x and y 
  #coordinates paired together
  
  #CREATE MATRIX OF DISTANCES BETWEEN PATCHES
  distance.matrix<-dist(coordinates[,1:2], method="euclidean", diag=TRUE, upper=TRUE)
  #calculate distances between patches based on euclidean distances. 
  #obvs dii=0 and dij=dji, which is indicated by diag=TRUE and upper=TRUE
  distance.matrix<-as.matrix(distance.matrix) #set d as a matrix
  
  #LOOP TO PREVENT PATCH OVERLAP
  for (i in 1:n.patches){ #for every patch
    for (j in 1:n.patches){ #and every other patch
      while (j!=i & distance.matrix[i,j] < (radii[i]+radii[j])) { 
        #while the distance between the 2 patches is less than the sum of their radii
        if (landscape.type == "2D"){ #if it's a 2D landscape
          x.coord[i]<-runif(1, min=0, max=landscape.limit) 
          #pick a new x coordinate for that patch
        } else if (landscape.type == "linear") { #otherwise if it's linear
          p.x.coord<-0} #just leave this as 0
        y.coord[i]<-runif(1, min=0, max=landscape.limit) 
        #pick a new y coordinate for that patch
        coordinates[i,]<-c(x.coord[i], y.coord[i]) #update the coordinates 
        distance.matrix<-dist(coordinates[,1:2], method="euclidean", diag=TRUE, upper=TRUE) 
        #update the distances
        distance.matrix<-as.matrix(distance.matrix) #set d as a matrix
      }}}
  
  if (patch.distribution != "random"){ #If we want a non-random landscape (clustered or regular)...
    for (i in 1:clustering.iters){ #for as many interations as specified by the number of clustering
      #algorithm iterations
      
      #STEP 1: PICK A RANDOM PATCH UP FROM THE LANDSCAPE
      r<-sample(1:n.patches,1, replace=T) #pick a random number between 1 and n.patches
      r.coordinates<-coordinates[-r,] 
      #remove that patch from the landscape (pick it up to be potentially relocated)
      
      #STEP 2: CALCULATE THE DISTANCE AND CONNECTIVITY OF THE PATCH CHOSEN TO OTHER PATCHES
      d.r.to.j<-rep(NA, (n.patches-1))
      connectivity.r.to.j<-rep(NA, (n.patches-1))
      for (j in 1:(n.patches-1)) { #for every patch in the network calculate...
        d.r.to.j[j] <- sqrt((x.coord[r]-r.coordinates[j,1])^2+(y.coord[r]-r.coordinates[j,2])^2) 
        #the distance between point x and all patches
        connectivity.r.to.j[j]<-exp(-d.r.to.j[j])} 
      #the connectivity of point x and all others
      connectivity.r<-sum(connectivity.r.to.j) 
      #the connectivity of point x is the sum of it's connectivity to all patches
      
      #STEP 3: PICK A RANDOM POINT
      if (landscape.type == "2D"){ #if it's a 2D landscape
        p.x.coord<-runif(1, min=0, max=landscape.limit) 
        #pick a random x coordinate for the point between 0 and the landscape extent
      } else if (landscape.type == "linear") { #otherwise if it's linear
        p.x.coord<-0} #just leave this as 0
      p.y.coord<-runif(1, min=0, max=landscape.limit) 
      #pick a random y coordinate for the point between 0 and the landscape extent
      
      #STEP 4: CALCULATE THE DISTANCE AND CONNECTIVITY OF THE NEW POINT TO OTHER PATCHES
      d.x.to.j<-rep(NA, (n.patches-1))
      connectivity.x.to.j<-rep(NA, (n.patches-1))
      for (j in 1:(n.patches-1)) { #for every patch in the network calculate...
        d.x.to.j[j] <- sqrt((p.x.coord-r.coordinates[j,1])^2+(p.y.coord-r.coordinates[j,2])^2) 
        #the distance between point x and all patches
        connectivity.x.to.j[j]<-exp(-d.x.to.j[j])} 
      #the connectivity of point x and all others
      connectivity.x<-sum(connectivity.x.to.j) 
      #the connectivity of point x is the sum of it's connectivity to all patches
      
      #STEP 5: LOOP TO ENSURE NEW POINT WILL NOT RESULT IN PATCH OVERLAP
      for (g in 1:(n.patches-1)) { #for each patch
        while (g!=r & d.x.to.j[g]<(radii[r]+radii[g])){ 
          #while any patch other than the chosen patch r has a distance between it's 
          #radius and r's radius greater than the distance between that patches center 
          #and the proposed new location for r
          if (landscape.type == "2D"){ #if it's a 2D landscape
            p.x.coord<-runif(1, min=0, max=landscape.limit) 
            #pick a new x coordinate for the new point
          } else if (landscape.type == "linear") { #otherwise if it's linear
            p.x.coord<-0} #just leave this as 0
          p.y.coord<-runif(1, min=0, max=landscape.limit) 
          #pick a new y coordinate for the new point
          d.x.to.j<-rep(NA, (n.patches-1))
          for (j in 1:(n.patches-1)){ 
            #for every patch in the network update calculations for...
            d.x.to.j[j]<-sqrt((p.x.coord-r.coordinates[j,1])^2+(p.y.coord-r.coordinates[j,2])^2) 
            #the distance between point x and all patches
            connectivity.x.to.j[j]<-exp(-d.x.to.j[j])} 
          #the connectivity of point x and all others
          connectivity.x<-sum(connectivity.x.to.j) 
          #the connectivity of point x is the sum of it's connectivity to all patches
        }}
      
      #STEP 6: COMPARE THE CONNECTIVITY OF THE PATCH AT ITS ORIGINAL LOCATION TO ITS 
      #NEW LOCATION
      #to generate clustered landscape: increase the pr(r relocated) 
      #if connectivity.x > connectivity.r
      if (patch.distribution == "clustered"){ #if we wish to create a clustered landscape...
        if (connectivity.x > connectivity.r){ #if the new location has a higher connectivity...
          if (landscape.type == "2D"){ #if it's a 2D landscape
            x.coord[r]<-p.x.coord #set the x coordinate of the patch to be the new location's
          } else if (landscape.type == "linear") { #otherwise if it's linear
            p.x.coord<-0} #just leave this as 0
          y.coord[r]<-p.y.coord #set the y coordinate of the patch to be the new location's
          coordinates[r,]<-c(x.coord[r],y.coord[r]) 
          #enter the patch's new coordinates into the landscape
        }}
      #to generate regular landscape: decrease the pr(r relocated) 
      #if connectivity.x > connectivity.r
      if (patch.distribution == "regular") { #if we wish to create a regular landscape...
        if(connectivity.x < connectivity.r) { 
          #if the new location has a lower connectivity...
          if (landscape.type =="2D"){ # if it's a 2D landscape
            x.coord[r]<-p.x.coord #set the x coordinate of the patch to be the new location's
          } else if (landscape.type == "linear") { #otherwise if it's linear
            p.x.coord<-0} #just leave this as 0
          y.coord[r]<-p.y.coord #set the y coordinate of the patch to be the new location's
          coordinates[r,]<-c(x.coord[r],y.coord[r]) 
          #enter the patch's new coordinates into the landscape
        }}
    }
  }
  } else if (patch.distribution=="uniform") { #if it's a uniform (grid in 2D) patch distribution
    if (landscape.type == "2D"){ #if it's 2D
      x.coord<-landscape.limit/n.patches*patch.ID-(landscape.limit/n.patches)/2 
      #set the x.coords to be evenly spaced throughout the landscape limit
    } else if (landscape.type=="linear"){ #if it's linear
      x.coord<-rep(0, n.patches) #arrange the patches along the y axis (x=0)
    }
    y.coord<-landscape.limit/n.patches*patch.ID-(landscape.limit/n.patches)/2
    #and set the y.coords to be evenly spaced within the landscape limit
    coordinates<-data.frame(x.coord,y.coord)
  }
  
  #STEP 7: UPDATE MATRIX OF DISTANCES BETWEEN PATCHES
  distance.matrix<-dist(coordinates[,1:2], method="euclidean", diag=TRUE, upper=TRUE)
  #calculate distances between patches based on euclidean distances. 
  #obvs dii=0 and dij=dji, which is indicated by diag=TRUE and upper=TRUE
  distance.matrix<-as.matrix(distance.matrix) #set d as a matrix
  
  landscape<-data.frame(patch.ID, areas, coordinates, distance.matrix)
  return(landscape)}
#END OF CREATE LANDSCAPE FUNCTION
#######################################################################################
#CHECKS:
#library("ggplot2")
#landscape<-create.landscape(n.patches=4, landscape.type="2D", landscape.limit=20, 
#                 patch.distribution="clustered", areas.distribution="random", areas.limit=2, 
#                 clustering.iters=10)
#head(landscape)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))

#landscape<-create.landscape(n.patches=4, landscape.type="linear", landscape.limit=20, 
#                 patch.distribution="regular", areas.distribution="uniform", areas.limit=2, 
#                 clustering.iters=10)
#head(landscape)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))

#landscape<-create.landscape(n.patches=4, landscape.type="linear", landscape.limit=20, 
#                            patch.distribution="uniform", areas.distribution="uniform", 
#                            areas.limit=2, clustering.iters=10)
#head(landscape)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))

