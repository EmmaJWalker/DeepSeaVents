##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function creates unique names for each element of a matrix.

# INPUTS:
# mat: a matrix
# mat.name: a character string providing a name for the matrix (each element 
#           will be named by this followed by it's row and column)

# OUTPUTS: 
# m.names: a matrix of character strings providing names for each element of 
#          the matrix where each element will be named by this followed by it's 
#          row and column.

# REQUIRED PACKAGES:
# none

# REQUIRED FUNCTIONS: 
# none

####################################################################################
# FUNCTION CODE:
####################################################################################
name_Mij<-function(mat, mat.name, MorV){
  m.names<-matrix(rep(NA,nrow(mat)*ncol(mat)),nrow(mat),ncol(mat))
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      m.names[i,j]<-paste0(mat.name,i,j)
    }
  }
  return(m.names)
}
####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#mat<-diag(2)
#mat.name<-"I"
#name_Mij(mat=mat, mat.name=mat.name)
##########################################################################################################


