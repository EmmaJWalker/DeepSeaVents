##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function extracts the distance matrix from a landscape 

# INPUTS:
#landscape: a dataframe containing patch.ID, p.sizes, coordinates, distance.matrix
#            -patch.ID is a unique integer ID number for each patch
#            -p.sizes provides the size of each patch
#            -y.coord provides vertical (e.g. latitudinal) patch locations on a "map"
#            -x.coord provides horizontal (e.g. longitudinal) patch locations on a "map" (if 2D landscape)
#            -subsequent columns contain interpatch distances forming a distance (or adjacency) matrix

# OUTPUTS:
# dist.mat: a matrix of interpatch distances within the input landscape

# REQUIRED PACKAGES:
# none

####################################################################################
# FUNCTION CODE:
####################################################################################

get_distmat<-function(landscape){
  dist.mat<-landscape[,-(1:4)]
  return(dist.mat)
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("ggplot2")
#landscape<-make_landscape(n.patches=4, landscape.type="plane", landscape.size=20, 
#                 p.dist="clustered", p.size.dist="random", max.p.size=2, 
#                 clust.its=10)
#head(landscape)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = p.sizes)))
#dist.mat<-get_distmat(landscape)
#dist.mat
####################################################################################