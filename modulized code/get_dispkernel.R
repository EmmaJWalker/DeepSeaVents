##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function calculates a dispersal kernel as it might be specified for a species' metapopulation

# INPUTS:
# dist.mat: a matrix of interpatch distances 
# alpha: a positive numeric value indicating 1/(avg. dispersal distance of a species)
# gamma: a positive numeric value <1 indicating the strength to which upstream (where upstream are those with j>i
#        as indexed by patchID when "landscape" is ordered by patchID) dispersal may be limited (0 no limitation)
# self.rec: a positive numeric value indicating the extent to which disperses from a given patch may resettle on 
#        the same patch (0 makes resettlement impossible)

# OUTPUTS:
# disp.kernel: a dataframe of values weighted relative to a species ability to move from one patch 
#              location to another (i.e. effective adjacency of patches)

# REQUIRED PACKAGES:
# none

####################################################################################
# FUNCTION CODE:
####################################################################################
get_dispkernel<-function(dist.mat, alpha, gamma, self.rec){
  n.patches<-nrow(dist.mat)
  k<-matrix(rep(0,n.patches*n.patches),n.patches,n.patches)
  
  for (i in 1:n.patches){ #for each patch
    for(j in 1:n.patches){ #to each other patch
      if(i<j){ #if patch j is "upstream" of patch i
        k[i,j]<-(1-gamma)*exp(-alpha*dist.mat[i,j])} #dispersal decreases according to the strength
      #to which upstream dispersal is limited, and further decreases exponentially with 
      #1/(avg.dispersal distance) (alpha) and the distance between the patches
      else if (i>j) { #if patch j is "downstream" of patch i
        k[i,j]<-exp(-alpha*dist.mat[i,j])} #dispersal decreases exponentially with 
      #1/(avg. dispersal distance) (alpha) and the distance between the patches
      else if (i==j){ #otherwise if the other patch is itself and self recruitment
        #is allowed
        k[i,j]<-1} #dispersal to itself occurs
      else{ #otherwise,
        k[i,j]<-0 #no dispersal between the two patches or itself occurs
      }
    }
  }
  return(k) #provide the appropriate dispersal kernel
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
#disp.kernel<-get_dispkernel(dist.mat, alpha=0.1, gamma=0.5, self.rec=1)
#disp.kernel
####################################################################################
