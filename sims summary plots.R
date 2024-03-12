rm(list=ls()) #clear the workspace
plot.data<-read.csv("plot_data_test1.csv")
head(plot.data)

library(ggplot2)
library(wesanderson)
library(ghibli)
wes_palettes$Darjeeling1
ghibli_palettes


plot.data$r.off<-as.factor(log10(plot.data$r.off))
plot.data$r.on<-as.factor(log10(plot.data$r.on))
plot.data$alpha<-factor(plot.data$alpha, 
                        levels=c("0.01", "0.1", "10"),
                        labels=c("Global Dispersal", "Stepping Stone", "1/10th Stepping Stone"))
plot.data$percent.occupied<-plot.data$avg.sizes/10*100




# grouped boxplot avg. metapop size
ggplot(plot.data) + 
  geom_boxplot(aes(x=r.on, y=exp.QED.size/10*100, fill=r.off), col="red", alpha=0.75) + 
  geom_boxplot(aes(x=r.on, y=delta, fill=r.off), col="blue", alpha=0.75) + 
  geom_boxplot(aes(x=r.on, y=percent.occupied, fill=r.off), col="black", alpha=0.75) + 
  theme_classic() + scale_fill_manual(values=ghibli_palettes$MononokeMedium[-1]) +
  labs(title = "Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes",
       subtitle = "(measured in % habitat occupied over 10,000 simulated timesteps)",
       y = "Average % Habitat Occupied", x = "Log Rate of Recovery") + 
  guides(fill=guide_legend("Log Rate of Disturbance")) +
  facet_grid(~ alpha)

# grouped boxplot avg. number of habitat patches
ggplot(plot.data) + 
  geom_boxplot(aes(x=r.on, y=exp.n.on, fill=r.off), col="red", alpha=0.75) +
  geom_boxplot(aes(x=r.on, y=avg.n.ons, fill=r.off), col="black", alpha=0.9) + 
  theme_classic() + scale_fill_manual(values=ghibli_palettes$MarnieMedium2[-2]) +
  labs(title = "Expected versus Average Amount of Habitat Across Disturbance Regimes",
       subtitle = "(measured in number of habitatable patches over 10,000 simulated timesteps)",
       y = "Average Number of Habitatable Patches", x = "Log Rate of Recovery") + 
  guides(fill=guide_legend("Log Rate of Disturbance"))

######################################## some fun with colours ######################
#ggplot(plot.data[alpha==0.01,]) + 
#  geom_boxplot(aes(x=r.on, y=avg.sizes, fill=r.off), col="white") + 
#  scale_color_manual(values=ghibli_palettes$MonokeMedium) +
#  theme_dark() + scale_fill_manual(values=ghibli_palettes$MononokeLight) #just testing

#ggplot(plot.data[alpha==0.1,]) + 
#  geom_boxplot(aes(x=r.on, y=avg.sizes, fill=r.off), col="black", alpha=0.9) + 
#  theme_classic() + scale_fill_manual(values=ghibli_palettes$PonyoMedium[-1]) #also nice

#ggplot(plot.data[alpha==0.01,]) + 
#  geom_boxplot(aes(x=r.on, y=avg.n.ons, fill=r.off), col="black", alpha=0.75) + 
#  theme_classic() + scale_fill_manual(values=ghibli_palettes$YesterdayDark[-2])

