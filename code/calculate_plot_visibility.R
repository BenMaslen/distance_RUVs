library(dplyr)
library(ggplot2)


dist_prob_fun <- function(x){
  dnorm(x, mean = 0, sd = imp_sd, log = FALSE)/dnorm(0, mean = 0, sd = imp_sd, log = FALSE)
}

#we want 20m vis

#one plot

norm_x = 0:30

vis_perc = 0.05
dist = 3

imp_sd = dist/sqrt(-2*log(vis_perc))

norm_y = dist_prob_fun(norm_x)
imp_func_dat <- data.frame(norm_x,norm_y)

imp_func_dat[21,]

ggplot(imp_func_dat,aes(x=norm_x,y=norm_y)) + geom_line() +
  xlab("Distance") + ylab("Probability of detection") +
  theme_classic()

#all lines



norm_x = seq(0,30,by=0.001)

vis_perc = 0.05

detection <- vector()
distance <- vector()

for (i in 5:20){
  dist = i
  imp_sd = dist/sqrt(-2*log(vis_perc))
  detection = c(detection,dist_prob_fun(norm_x))
  distance = c(distance,rep(i,length(norm_x)))
}

detection_dat = data.frame(detection,distance,norm_x=rep(norm_x,length(5:20)))



ggplot(detection_dat,aes(x=norm_x,y=detection,group=distance,colour=distance)) + geom_line() +
  xlab("distance (m)") + ylab("probability of detection") +
  theme_classic() + scale_color_continuous(name="visibility (m)")
  
ggsave("plots/figure_3.png",height=1.1*5.12*0.8,width=1.1*5.79*(0.98/0.62)*0.8)



