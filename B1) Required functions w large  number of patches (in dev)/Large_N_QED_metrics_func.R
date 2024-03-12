#USE ONLY A SUBSAMPLE OF THE POSSIBLE CONFIGURATIONS:
configs.sample<-sample((2:(2^n.patches)), size=100, replace=F)
int.To.Bin(configs.sample)