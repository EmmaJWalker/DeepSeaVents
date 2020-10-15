#Figure 6 compares the calculated metapopulation size at q.e.d. (red
#line) to a histogram of mean population sizes from a number of
#simulations of the combined model, with mean value given by the black
#line. This set of parameters lies within the region for which the
#metapopulation size at q.e.d. gives a good approximation of the actual
#time-averaged metapopulation size.


#This script uses the quasi-equilibrium distribution to obtain the
#metapopulation size at quasi-equilibrium (by weighting the q.e.d. with
#the equilibrium metapopulation size for each configuration), and then
#compares that value to the actual time-averaged metapopulation size
#for a given simulation of the CTMC+SRLM, obtained by integrating the
#simulated metapopulation size over time and dividing by the time interval.
rm(list=ls()) #clear the workspace
library("pracma")
library("deSolve")
setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned")
source("create_landscape_func.r")
source("SRLM_ODE_func.r")
source("at_equilibrium_rootfunc.r")
source("dispersal_kernel_func.r")

#Initializing Parameters:
###################################################################################################
total.t<-5000 #total time to run the model (optional; could run to absorption)
n.patches<-6 #number of sites
A<-runif(n.patches, min=0, max=4) #Instead of scaling by a certain area we define the interval
locs<-create.landscape(6, "linear", landscape.limit=10)
e.rate<-0.3
c.rate<-0.5
delta<-e.rate/c.rate
alpha<-0.1 #inverse of the mean dispersal distances (negetive exponential)
#network structure parameters:
gamma<-0
epsilon<-1
self.rec<-1
#rates of vents turning off to on and on to off respectively
r1<-0.02
r2<-0.01
config<-rbinom(n.patches, 1, 0.5)
####################################################################################################

#calculations to obtain all our parameters
####################################################################################################
#plot the habitat network and calculate necessary preliminary matrices/values
x.coord<-locs[,1]
y.coord<-locs[,2]
#calculate E for each patch
E<-e.rate/A
# construct the distance matrix
dist.mat<-matrix(rep(0,n.patches*n.patches),n.patches,n.patches)
for (i in 1:n.patches){
  for (j in 1:n.patches){
    dist.mat[i,j]<-sqrt((x.coord[i]-x.coord[j])^2+(y.coord[i]-y.coord[j])^2)
  }
}
#Getting the f(dij) values to plug into the rhs
f.vals<-disp.kernel(dist.mat,alpha,gamma,epsilon,self.rec,n.patches)
f.vals<-as.vector(f.vals)
f.names<-matrix(rep(NA,n.patches*n.patches),n.patches,n.patches)
for (i in 1:n.patches){
  for (j in 1:n.patches){
    f.names[i,j]<-paste0("f.vals",i,j)
  }
}
f.names<-as.vector(f.names)
###################################################################################################

runs<-1000
avg.size<-rep(0,runs)

#options1 = odeset('RelTol',1e-8,'AbsTol',1e-8); ????
#timer vectors to figure out how long things take
vec.times<-rep(0,runs)
sim.times1<-rep(0,runs)
sim.times2<-rep(0,runs)
sim.times3<-rep(0,runs)
for (run.val in 1:runs){
  
  #creating a list of the parameters we calculated to be passed to our ODE function
  ###################################################################################################
  parameter.names<-c("n.patches", "c.rate", "self.rec",
                     rep(paste0("E",1:n.patches)),
                     rep(paste0("x.coord",1:n.patches)),
                     rep(paste0("y.coord",1:n.patches)),
                     rep(paste0("A",1:n.patches)),
                     rep(paste0("config",1:n.patches))
                     ,f.names
  )
  parameter.values<-c(n.patches,c.rate,self.rec,E,x.coord,y.coord,A,config
                      ,f.vals
  )
  parameters<-list(c(setNames(parameter.values, parameter.names)))
  ###################################################################################################
  
  #setting our initial values
  ###################################################################################################
  new.ICs<-0.5*config
  IC.names<-rep(paste0("p",1:n.patches))
  IC.values<-new.ICs
  p<-setNames(IC.values,IC.names)
  tau<-4000 #just for checking
  ###################################################################################################
  
  for (k in 1:length(config)){
    
  }
}
output<-ode(SRLM.ODE, parameters,times=seq(0,round(tau,0),1),y=p)
                       N);



for k = 1:length(init_config)
vec(length(init_config) - k + 1) = str2double(init_config(k));
end
vec_times(run_val) = toc;

% this section starts the simulation process
t = 0;
% setting the very beginning initial condition to be half-occupied
sites at
% those sites randomly picked to be turned on above
new_ICs = 0.5*vec;
while (sum(vec) ~= 0)
  %while ((t < T) && (sum(vec) ~= 0))
    %first, figure out how long the CTMC stays in the current
state, and
%run the SRLM for that amount of time
rate_vec = zeros(N,1);
tic;
for i = 1:N
if (vec(i) == 0)
  rate_vec(i) = r1;
else
  rate_vec(i) = r2;
end
end
run_sum = sum(rate_vec);
tau = exprnd(1/run_sum);
left_out = find(~vec);
sim_times1(run_val) = toc;
tic;
[tee, p] = ode45(@(t,p)
                 SRLM_ode_TD_selfrec(t,p,N,c_rate,E,f_vals,...
                                     areas,left_out,self_rec),[0 tau],new_ICs,options1);
sim_times2(run_val) = toc;
tic;
tee = tee + t;
t_vals = [t_vals; tee];
p_vals = [p_vals; p];
sim_times3(run_val) = toc;
unif_val = rand;
comp_val = cumsum(rate_vec/run_sum);
for j = 1:N
if(unif_val <= comp_val(j))
  pos = j;
break

end
end
t = t + tau;
if (vec(j) == 0)
  vec(j) = 1;
elseif (vec(j) == 1)
vec(j) = 0;
end
new_ICs = zeros(N,1);
new_ICs(find(vec)) = p_vals(end,find(vec));
end
p_vals_weighted = zeros(length(p_vals),1);
for i = 1:length(p_vals)
p_vals_weighted(i) = p_vals(i,:)*areas/(sum(areas));
end
%plot(t_vals,p_vals_weighted)
%hold on;
avg_size(run_val) = trapz(t_vals,p_vals_weighted)/(t_vals(end));
%avg_size(run_val)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now calculate the metapopulation size at quasi-equilibrium and
compare
% go through each configuration and figure out whether it is stable,
and
% figure out the equilibrum metapopulation size
M_mat = zeros(N,N);
for i = 1:N
for j = 1:N
M_mat(i,j) =
  f_vals(i,j)*areas(i)*areas(j);
end
end
lambda_M = zeros(2^N,1);
metapop_size = zeros(2^N,1);
for i = 0:(2^N-1)
bin_val = dec2bin(i);
vec = zeros(N,1);
for k = 1:length(bin_val)
vec(length(bin_val) - k + 1) = str2double(bin_val(k));

end
left_out = find(~vec);
new_ICs = zeros(N,1);
new_ICs(find(vec)) = 0.5;
if(i == 0)
  metapop_size(i+1) = 0;
else
  % when i = 0, it's the config where everything is turned off
so
% pulling out the matrix elements doesn't work
M_mat_config = M_mat(find(vec),find(vec));
[A B C] = eig(M_mat_config);
max_index = find(max(B) == max(max(B)));
lambda_M(i+1) = max(max(B));
if(lambda_M(i+1) < delta)
  metapop_size(i+1) = 0;
else
  % run the ode until it hits equilibrium, to get the
equilibrium
% metapopulation size for the configuration
odefun = @(t,p)
SRLM_ode_TD_selfrec(t,p,N,c_rate,E,f_vals,areas,left_out,self_rec);
efun = @(t,p)
eventfun(t,p,N,c_rate,E,f_vals,areas,left_out,self_rec);
options = odeset('Events',efun);
[t,p] = ode15s(odefun,[0 3000], new_ICs, options);
end_p = p(end,:);
metapop_size(i+1) =
  dot(areas,end_p)/sum(areas);
end
end
end
%%
  %[dom_left, T_absorp] = quasi_eq_dist(N,r1,r2);
[dom_left, T_absorp] = qed_calculator(N,r1,r2);
weighted_metapop_size = dot(dom_left,metapop_size(2:end));
histogram(avg_size(avg_size >= 0),25,'Normalization','probability')
line([weighted_metapop_size, weighted_metapop_size], [0,1])
% % extra code to play around with the run times
%
% y1 = vec_times;
% y2 = y1 + sim_times1;
% y3 = y2 + sim_times2;
% y4 = y3 + sim_times3;