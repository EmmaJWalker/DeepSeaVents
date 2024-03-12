##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# Assuming you have an indexed set of jobs (numbered 1-J) to be run over only N computational cores, 
# this function splits that set into an indexed subset of those jobs to be run on a core. 
# e.g. if I have 10 jobs numbered 1-10 to run in parallel and only 3 computational cores available
# this function will give jobs 1, 2, 3, and 4 to core 1, jobs 5, 6, and 7 to core 2, and jobs 8, 9, 10
# to core 3.
# N is the maximum number of core you would like the jobs split over. If the number of jobs (J) is less
# than N, unnecessary cores will remain unused.
# if I had 2 jobs, core gets 1 job, core 2 gets job 2, and core 3 is unused.

# INPUTS:
# N: the maximum number of cores to divide the jobs over
# J: the number of jobs to divide over them
# i: the core index (i.e. index indicating which core to use)

# OUTPUTS:
# jobs.for.i: a vector providing an index of the jobs to be run on core i

# REQUIRED PACKAGES:
# none

# REQUIRED FUNCTIONS:
# none

####################################################################################
# FUNCTION CODE:
####################################################################################
split_jobs<-function(J, N, i){
  if (J<N){ #if the number of jobs to be given is less than or equal to the number of cores
    #give out 1 job to the first J cores
    jobs<-1:J # thus, for each job
    if (i<=J){ # if i<=J
      jobs.for.i<-jobs[i] #give it one of those jobs
    } else {# else do nothing
    jobs.for.i<-NA
      }
  } else { #otherwise if there are more jobs than cores
    v<-floor((J)/N) #number of jobs that can be evenly divided to each core
    remainder<-J-(N*v) #set of remaining jobs
    jobs.for.i<-c((((i-1)*v)+1):(i*v)) #vector of jobs as can be evenly assigned to each core
    if (remainder>0){ # if there are any remainling jobs that couldn't be evenly distributed to each core
      if (i<=remainder){ #and if i is less than the remainder
        jobs.for.i<-c(jobs.for.i,(i*v+1)) #add another job to i
        jobs.for.i<-jobs.for.i+(i-1) #shift the job indices appropiately to account for the remaining jobs
        # that were added
      } else {
        jobs.for.i<-jobs.for.i+remainder
      }
    }
  }
  return(jobs.for.i)
}


####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#J<-10
#N<-3
#for(i in 1:N){
#  jobs.for.i<-split_jobs(J=J,N=N, i=i)
#  print(jobs.for.i)
#  print(i)
#}
#J<-2
#N<-3
#for(i in 1:N){
#  jobs.for.i<-split_jobs(J=J,N=N, i=i)
#  print(jobs.for.i)
#  print(i)
#}
####################################################################################
