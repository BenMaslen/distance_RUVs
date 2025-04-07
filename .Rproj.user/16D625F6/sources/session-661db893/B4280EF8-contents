library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

fig_6_dat <- read.csv("data/Fig_6_dat.csv")



ggplot(fig_6_dat,aes(y=D+100,x=MaxN+0.1,colour=Species,shape=Site))+geom_point() + theme_classic() +
  scale_y_log10() + scale_x_log10()
ggplot(fig_6_dat,aes(y=D+100,x=MeanCount+0.1,colour=Species,shape=Site))+geom_point() + theme_classic() +
  scale_y_log10() + scale_x_log10() 


fig_6a = ggplot(fig_6_dat,aes(y=D,x=MaxN,colour=Species,shape=Site))+geom_point() + theme_classic() + ylab("Distance sampling density")
fig_6b = ggplot(fig_6_dat,aes(y=D,x=MeanCount,colour=Species,shape=Site))+geom_point() + 
  theme_classic() +
  theme(#axis.line.y=element_blank(),
       #axis.text.y=element_blank(),axis.ticks=element_blank(),
       axis.title.y=element_blank())



ggarrange(fig_6a,fig_6b,ncol=2,common.legend = T)
