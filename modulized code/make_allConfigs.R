# This function constructs a sparse matrix where each row represents each unique configuration a landscape with
# a given number of patches could take on, where a 1 or 0 in each column indicates whether a patch identified by
# that column being either present or absent within a given landscape of those patches. 

# The row order of configurations corresponds to their row and column order in the generator matrix (G_mat)
# and likewise, the construction of habitat configurations takes advantage of block structure by which
# the generator matrix grows to save on computation and memory - see the Configuration IDs document for details.

# While this function is extremely computationally and memory efficient, the task completed by this function becomes 
# impossible for landscapes with large numbers of patches (n.patches) because the number of configurations possible 
# grows to 2^n.patches and thus becomes huge very quickly! 

# INPUTS:
# n.patches: the number of patches in a landscape

# OUTPUTS:
# a binary sparse matrix where each row represents each unique configuration a landscape with
# a given number of patches could take on, where a 1 or 0 in each column indicates whether a patch identified by
# that column being either present or absent within a given landscape of those patches.

# REQUIRED PACKAGES:
# Matrix: to get sparse matrices and vectors

# REQUIRED FUNCTIONS:
# none
##########################################################################################################
make_allConfigs<-function(n.patches){
  n<-n.patches
  all.configs <- Matrix(0, nrow=2^n, ncol=n, sparse=TRUE) #sparse matrix version
  #all.configs<-matrix(rep(NA, (2^n)*n), nrow=2^n, ncol=n) #no sparse matrix version
  for (i in 1:n){
    #if (i==1){BB<-c(0,1) #no sparse vector version
    if (i==1){BB<-as(c(0,1), "sparseVector") #sparse vector version
    }else{BB<-c(BB[1:(2^(i-2))],BB,BB[((2^(i-2))+1):(2^(i-1))])}
    all.configs[,(n-i+1)]<-t(rep(BB,2^(n-i)))}
  return(all.configs)
}
##########################################################################################################
# QUICK CHECK: (uncomment and run to check)
##########################################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("Matrix") #to get sparse matrices
#configs<-make_allConfigs(n.patches=15) # 10 takes about a seccond, 15 several secconds, ...
#configs
##########################################################################################################

