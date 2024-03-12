##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function creates a landscape of locations at which local habitat patches may be formed ("turn on") 
# or lost ("turn off").

# There are parameters to control:
# 1. how many patches may be present (i.e. the number of locations at which individual patches are formed 
#    and lost) 
#      - this is simply given by "n.patches"
# 2. how the sizes (i.e. p.sizes or carrying capacities) of these patches vary within the landscape
#      - by setting how patch sizes are distributed within the landscape ("lognormal"-ly, "uniform"-ly or 
#        "random"-ly) for "p.size.dist"
#        and or/a maximum size ("max.p.size") to avoid patches exceeding the desired size of the entire 
#        landscape
# 3. how big to make the entire landscape 
#     - use "landscape.size" to indicate the max x and y coordinates for where patches can be 
# 3. whether patches are arranged along a 1D "line" or within a 2D "plane"
#     - specify "line" or "plane" for "landscape.type"
# 4. whether the patches tend to be "clustered" more together, versus "random"-ly or more "uniformal"-ly 
#    distributed throughout the landscape
#     - specify "clustered", "random", or "uniform" for "p.dist"
#     - adjusting "clust.its" to be a higher integer will increase the clustering


# INPUTS:
# n.patches: integer equal to the total number of patches in the system
# landscape.type: "line" or "plane" indicating whether to have a 1D or 2D patch arrangement
# landscape.size: numerical value indicating max y (and x coordinates if 2D)
# p.dist: "clustered", "random", "uniform"
# clust.its: integer equal to the number of times to iterate the algorithm increasing clustering or uniformity for
#            the distribution of patches in the landscape
# p.size.dist: "lognormal", "random" or "uniform" indicating how patch sizes are to be distributed
# max.p.size: numerical value indicating largest patch size

# OUTPUTS:
# landscape: a dataframe containing patch.ID, p.sizes, coordinates, distance.matrix
#            -patch.ID is a unique integer ID number for each patch
#            -p.sizes provides the size of each patch
#            -y.coord provides vertical (e.g. latitudinal) patch locations on a "map"
#            -x.coord provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
#            -subsequent columns contain interpatch distances forming a distance (or adjacency) matrix

# REQUIRED PACKAGES:
# none

####################################################################################
# FUNCTION CODE:
####################################################################################
make_landscape<-function(n.patches, landscape.type, landscape.size, p.dist, 
                             p.size.dist, max.p.size, clust.its){
    n.patches<-n.patches
    patch.ID<-c(1:n.patches)
    #assign ID numbers to each patch (This is so that we can track which patches are
    #removed throughout destruction of patches from the landscape)
    if (p.size.dist == "lognormal") { #if p.sizes are to be distributed lognormally
      p.sizes<-rlnorm(n.patches, meanlog = 2, sdlog=1 )
      #patch p.sizes lognormally distributed with a mean of 2 and standard deviation of 1
      #this was chosen to best mimic most natural landscapes with lognormally distributed p.sizes
      #but also set so that we don't get one patch taking up almost all the landscape 
      #(acting like a mainland, that every other patch's dynamics depends very heavily on)
      #within a 100x100 unit landscape
    }
    if (p.size.dist == "uniform"){ #if p.sizes are to be distributed uniformly
      p.sizes<-rep(max.p.size, n.patches) #patch p.sizes will all be at max.p.size
    }
    if (p.size.dist == "random"){ #if p.sizes are to be distributed randomly
      p.sizes<-runif(n.patches, min=1, max=max.p.size) #patch p.sizes will be randomly selected from a 
      #uniform distribution between one and the area limit
    }
    
    radii<-sqrt(p.sizes/pi) #calculate radii of hypothetically circular patches
    
    if (p.dist!="uniform"){#if not arranging aranging pactches uniformly (grid in plane)
      
      if (landscape.type == "line"){ #if creating a line landscape
        x.coord<-rep(0, n.patches) #patches will be arranged vertically along the y axis (x = 0)
      } else if (landscape.type == "plane"){
        x.coord<-runif(n.patches, min=0, max=landscape.size) 
        #pick a random number for the x coordinate of each patch between 0 and the extent
        #of the landscape
      }
      y.coord<-runif(n.patches, min=0, max=landscape.size) 
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
            if (landscape.type == "plane"){ #if it's a plane landscape
              x.coord[i]<-runif(1, min=0, max=landscape.size) 
              #pick a new x coordinate for that patch
            } else if (landscape.type == "line") { #otherwise if it's line
              p.x.coord<-0} #just leave this as 0
            y.coord[i]<-runif(1, min=0, max=landscape.size) 
            #pick a new y coordinate for that patch
            coordinates[i,]<-c(x.coord[i], y.coord[i]) #update the coordinates 
            distance.matrix<-dist(coordinates[,1:2], method="euclidean", diag=TRUE, upper=TRUE) 
            #update the distances
            distance.matrix<-as.matrix(distance.matrix) #set d as a matrix
          }}}
      
      if (p.dist != "random"){ #If we want a non-random landscape (clustered or regular)...
        for (i in 1:clust.its){ #for as many interations as specified by the number of clustering
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
          if (landscape.type == "plane"){ #if it's a plane landscape
            p.x.coord<-runif(1, min=0, max=landscape.size) 
            #pick a random x coordinate for the point between 0 and the landscape extent
          } else if (landscape.type == "line") { #otherwise if it's line
            p.x.coord<-0} #just leave this as 0
          p.y.coord<-runif(1, min=0, max=landscape.size) 
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
              if (landscape.type == "plane"){ #if it's a plane landscape
                p.x.coord<-runif(1, min=0, max=landscape.size) 
                #pick a new x coordinate for the new point
              } else if (landscape.type == "line") { #otherwise if it's line
                p.x.coord<-0} #just leave this as 0
              p.y.coord<-runif(1, min=0, max=landscape.size) 
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
          if (p.dist == "more clustered"){ #if we wish to create a clustered landscape...
            if (connectivity.x > connectivity.r){ #if the new location has a higher connectivity...
              if (landscape.type == "plane"){ #if it's a plane landscape
                x.coord[r]<-p.x.coord #set the x coordinate of the patch to be the new location's
              } else if (landscape.type == "line") { #otherwise if it's line
                p.x.coord<-0} #just leave this as 0
              y.coord[r]<-p.y.coord #set the y coordinate of the patch to be the new location's
              coordinates[r,]<-c(x.coord[r],y.coord[r]) 
              #enter the patch's new coordinates into the landscape
            }}
          #to generate regular landscape: decrease the pr(r relocated) 
          #if connectivity.x > connectivity.r
          if (p.dist == "regular") { #if we wish to create a regular landscape...
            if(connectivity.x < connectivity.r) { 
              #if the new location has a lower connectivity...
              if (landscape.type =="plane"){ # if it's a plane landscape
                x.coord[r]<-p.x.coord #set the x coordinate of the patch to be the new location's
              } else if (landscape.type == "line") { #otherwise if it's line
                p.x.coord<-0} #just leave this as 0
              y.coord[r]<-p.y.coord #set the y coordinate of the patch to be the new location's
              coordinates[r,]<-c(x.coord[r],y.coord[r]) 
              #enter the patch's new coordinates into the landscape
            }}
        }
      }
    } else if (p.dist=="uniform") { #if it's a uniform (grid in plane) patch distribution
      if (landscape.type == "plane"){ #if it's plane
        x.coord<-landscape.size/n.patches*patch.ID-(landscape.size/n.patches)/2 
        #set the x.coords to be evenly spaced throughout the landscape limit
      } else if (landscape.type=="line"){ #if it's line
        x.coord<-rep(0, n.patches) #arrange the patches along the y axis (x=0)
      }
      y.coord<-landscape.size/n.patches*patch.ID-(landscape.size/n.patches)/2
      #and set the y.coords to be evenly spaced within the landscape limit
      coordinates<-data.frame(x.coord,y.coord)
    }
    
    #STEP 7: UPDATE MATRIX OF DISTANCES BETWEEN PATCHES
    distance.matrix<-dist(coordinates[,1:2], method="euclidean", diag=TRUE, upper=TRUE)
    #calculate distances between patches based on euclidean distances. 
    #obvs dii=0 and dij=dji, which is indicated by diag=TRUE and upper=TRUE
    distance.matrix<-as.matrix(distance.matrix) #set d as a matrix
    
    landscape<-data.frame(patch.ID, p.sizes, coordinates, distance.matrix)
    return(landscape)
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("ggplot2")
#landscape<-make_landscape(n.patches=4, landscape.type="plane", landscape.size=20, 
#                 p.dist="more clustered", p.size.dist="random", max.p.size=2, 
#                 clust.its=10)
#head(landscape)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = p.sizes)))
####################################################################################




