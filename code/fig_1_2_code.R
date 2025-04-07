#Make figures 1 and 2

library(ggplot2)
library(latex2exp)
library(RColorBrewer)
library(dplyr)

#read in data

wa_dat = read.csv("results/fig_1_2_wa_dat.csv")
sim_dat = read.csv("results/sim_fig1_dat.csv")

#figure 2

#no reference point
ggplot(wa_dat,aes(plane_dist,height,colour=plane_dist>1.96757)) + 
  geom_point(alpha=0.6) + theme_classic() + ylim(0,0.8) + scale_x_continuous(limits=c(0,5),expand=c(0,0)) +
  geom_hline(yintercept = 0.78) + geom_vline(xintercept = 1.96757) + 
  geom_abline(intercept = 0.35,slope=0.2185765,linetype='dotted') + geom_hline(yintercept = 0.35,linetype='dotted') + 
  theme(legend.position="none") + xlab("distance (m)") + ylab("height (m)") + 
  geom_text(x = 0.27, y = 0.38, 
            label =  TeX("$\\hat{\\gamma}$"),size=3.5,colour="black")+
  geom_text(x = 1.9, y = 0, 
            label =  TeX("$\\hat{d}_{\\gamma}$"),size=3.5,colour="black")+
  geom_text(x = 0.17, y = 0.75, 
            label =  TeX("$\\upsilon\\hat{h}_{max}$"),size=3.5,colour="black")+
  geom_text(x = 0.14, y = 0.33, 
            label =  TeX("$h_{c}$"),size=3.5,colour="black")+
  scale_color_manual(values = c("#969696", "#252525"))

#with reference point
fig2 = ggplot(wa_dat %>% filter(Family_Genus_Species %in% c("Labridae:Choerodon:cyanodus")),aes(plane_dist,height,colour=plane_dist>1.762585)) + 
  geom_point(alpha=0.5) + theme_classic() + ylim(0,1.5) + scale_x_continuous(limits=c(0,5),expand=c(0,0)) +
  geom_hline(yintercept = 0.78) + geom_vline(xintercept = 1.762585) + 
  geom_abline(intercept = 0.35,slope=0.2439964,linetype='dotted') + 
  geom_abline(intercept = 0.35,slope=-0.3639964,linetype='dotted') + 
  geom_hline(yintercept = 0.35,linetype='dotted') + 
  theme(legend.position="none") + xlab("distance (m)") + ylab("height (m)") + 
  geom_text(x = 0.4, y = 0.4, 
            label =  TeX("$\\hat{\\gamma}_{+,c}$"),size=4,colour="black")+
  geom_text(x = 0.4, y = 0.3, 
            label =  TeX("$\\hat{\\gamma}_{-,c}$"),size=4,colour="black")+
  geom_text(x = 1.762585 + 0.16, y = 0, 
            label =  TeX("$\\hat{d}_{\\gamma_{+,c},j}$"),size=4,colour="black")+
  geom_text(x = 0.35/0.3639964 + 0.18, y = 0, 
            label =  TeX("$\\hat{d}_{\\gamma_{-,c}}$"),size=4,colour="black")+
  geom_text(x = 0.17, y = 0.72, 
            label =  TeX("$\\upsilon\\hat{h}_{max,j}$"),size=4,colour="black")+
  geom_text(x = 0.14, y = 0.31, 
            label =  TeX("$h_{c}$"),size=4,colour="black")+
  scale_color_manual(values = c("#969696", "#252525")) +
  geom_text(x = 4.209815+0.015, y = 1.37718, 
            label =  TeX("$X_U$"),size=4,colour="red")+
  geom_text(x = 0.5+0.015, y = 0.35-0.5*0.3639964, 
          label =  TeX("$X_L$"),size=4,colour="red") + 
geom_text(x = 4.209815+ 0.015 + 0.26, y = 1.37718 - 0.04, 
          label =  TeX("$(l_{\\gamma_{+,c}},h_{\\gamma_{+,c}})$"),size=3,colour="black")+
  geom_text(x = 0.5+0.015 + 0.26 , y = 0.35-0.5*0.3639964 - 0.04, 
            label =  TeX("$(l_{\\gamma_{-,c}},h_{\\gamma_{-,c}})$"),size=3,colour="black")



ggsave(fig2,"plots/figure_2.png",height=1.1*5.12*0.8,width=1.1*5.79*(0.98/0.62)*0.8)




 ggplot(wa_dat %>% filter(Family_Genus_Species %in% c("Labridae:Choerodon:cyanodus")),aes(plane_dist,height,colour=plane_dist>1.762585)) + 
  geom_point(alpha=0.5) + theme_classic() + ylim(0,1.5) + scale_x_continuous(limits=c(0,5),expand=c(0,0)) +
#  geom_hline(yintercept = 0.78) + geom_vline(xintercept = 1.762585) + 
  geom_abline(intercept = 0.35,slope=0.2439964,linetype='dotted') + 
  geom_abline(intercept = 0.35,slope=-0.3639964,linetype='dotted') + 
  geom_hline(yintercept = 0.35,linetype='dotted') + 
#  geom_hline(yintercept = 0.35,linetype='dotted') + 
  scale_color_manual(values = c("lightgray", "#252525")) +
  theme(legend.position="none") + xlab("distance (m)") + ylab("height (m)") 



ggsave("plots/figure_2_alt.png",height=1.1*5.12*0.8,width=1.1*5.79*(0.98/0.62)*0.8)


 ggplot(wa_dat %>% filter(Family_Genus_Species %in% c("Labridae:Choerodon:cyanodus")),aes(plane_dist-0.35,height,colour=plane_dist-0.35>1.762585)) + 
  geom_point(alpha=0.5) + theme_classic() + ylim(0,1.5) + scale_x_continuous(limits=c(0,5),expand=c(0,0)) +
#  geom_hline(yintercept = 0.78) + geom_vline(xintercept = 1.762585) + 
  geom_abline(intercept = 0.35,slope=0.2439964,linetype='dotted') + 
  geom_abline(intercept = 0.35,slope=-0.7139964,linetype='dotted') + 
  geom_hline(yintercept = 0.35,linetype='dotted') + 
#  geom_hline(yintercept = 0.35,linetype='dotted') + 
  scale_color_manual(values = c("lightgray", "#252525")) +
  theme(legend.position="none") + xlab("distance (m)") + ylab("height (m)") 



ggsave("plots/figure_2_alt2.png",height=1.1*5.12*0.8,width=1.1*5.79*(0.98/0.62)*0.8)


#figure 1

#prep and join data



wa_dat$type = "observed"
wa_dat$distance = wa_dat$plane_dist
sim_dat$type = "perfect detection"

fig_1_dat = rbind(wa_dat %>% dplyr::select(distance,type),sim_dat %>% dplyr::select(distance,type))

#fig_1_dat$type = factor(fig_1_dat$type,levels=c("perfect detection","observed"))

split_dat = fig_1_dat %>% group_by(type) %>% summarise(less2.5 = sum(distance<2.5),greater2.5 = sum(distance>2.5))

split_dat[2,2:3]/split_dat[1,2:3]


ggplot(fig_1_dat,aes(x=distance)) +
  geom_histogram(bins=20) + theme_classic() + xlab("distance (m)") + 
  facet_grid(type~.)

ggsave("plots/figure_1.png",height=1.1*5.12*0.8,width=1.1*5.79*(0.98/0.62)*0.8)



ggplot(fig_1_dat,aes(x=distance)) +
  geom_histogram(bins=5) + theme_classic() + xlab("distance (m)") + 
  facet_grid(type~.)

