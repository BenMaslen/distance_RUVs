#Distance sampling simulations

library(tidyverse)
library(MASS)
library(purrr)
library(Distance)
library(doRNG)

set.seed(12345)
overall_population_mean = 1952
#pop_disp = 4
N_vids = 13
R_max = 6
theta_diff = 51.99647
gamma_diff = 26.56505118
vid_length = 20*60
skipper = 30
N_stops = vid_length/skipper
z_max = 1
area = 100*100
N_sims = 1
imp_sd = 1.634156
speed = 10
schooling=1

abundance_estimates = rep(NA,N_sims)
sim_id = rep(NA,N_sims)


x_vids = runif(N_vids,min=R_max,max=sqrt(area)-R_max)
y_vids = runif(N_vids,min=R_max,max=sqrt(area)-R_max)
z_vids = 0
theta_min = runif(N_vids,min=-180,max=180-theta_diff)
theta_max = theta_min + theta_diff


vids_dat = data.frame(x_vids,y_vids,z_vids,theta_min,theta_max)


#N_pop_true = rnegbin(1,overall_population_mean,pop_disp)
N_pop_true = overall_population_mean/schooling


x = list()
y = list()
z = list()

for (k in (1:N_stops)){
  if(k==1){
    x[[k]] = runif(N_pop_true,min=0,max=sqrt(area))
    y[[k]] = runif(N_pop_true,min=0,max=sqrt(area))
    z[[k]] = runif(N_pop_true,min=0,max=z_max)
  }else{
    z[[k]] = z[[k-1]]
    x[[k]] = x[[k-1]] + rnorm(N_pop_true,0,speed/(sqrt(pi)/sqrt(2)))    #move in absolute value on average by speed m/m
    y[[k]] = y[[k-1]] + rnorm(N_pop_true,0,speed/(sqrt(pi)/sqrt(2)))    #move in absolute value on average by speed m/m
    
    for(l in (1:N_pop_true)){   #reflect points if we move outside the boundary
      if(x[[k]][l]<0){
        x[[k]][l] <- abs(x[[k]][l])
      }
      if(y[[k]][l]<0){
        y[[k]][l] <- abs(y[[k]][l])
      }
      
      if(x[[k]][l]>sqrt(area)){
        x[[k]][l] <- 2*sqrt(area)-x[[k]][l]
      }
      if(y[[k]][l]>sqrt(area)){
        y[[k]][l] <- 2*sqrt(area)-y[[k]][l]
      }
    }
  }
}

# Function to repeat each element in a vector 5 times
repeat_elements <- function(vec) {
  rep(vec, each = schooling)
}

# Apply the function to each vector in the list
x <- lapply(x, repeat_elements)
y <- lapply(y, repeat_elements)
z <- lapply(z, repeat_elements)





for (i in (1:N_vids)){
  for (k in (1:N_stops)){
    x_temp = x[[k]]
    y_temp = y[[k]]
    z_temp = z[[k]]
    
    pop_dat = data.frame(x=x_temp,y=y_temp,z=z_temp)
    
    d_temp_abs = sqrt((vids_dat$x_vids[i]-pop_dat$x)^2 + (vids_dat$y_vids[i]-pop_dat$y)^2 + (pop_dat$z)^2)
    d_temp_flat = sqrt((vids_dat$x_vids[i]-pop_dat$x)^2 + (vids_dat$y_vids[i]-pop_dat$y)^2)
    
    
    theta_temp = atan2((pop_dat$x - vids_dat$x_vids[i]),(pop_dat$y - vids_dat$y_vids[i]))*180/pi
    gamma_temp = atan(pop_dat$z/(sqrt((vids_dat$x_vids[i]-pop_dat$x)^2 + (vids_dat$y_vids[i]-pop_dat$y)^2)))*180/pi
    in_vid = d_temp_abs <= R_max & theta_temp>=theta_min[i] & theta_temp<=theta_max[i] & gamma_temp <= gamma_diff
    in_angle = theta_temp>=theta_min[i] & theta_temp<=theta_max[i] & gamma_temp <= gamma_diff
    
    if(i==1&k==1){
      dist_dat = data.frame(Sample.Label=i,step_id=k,distance=d_temp_flat,
                            distance_abs=d_temp_abs,Effort=N_stops,in_vid,
                            in_angle,Region.Label="Sims",Area=area,N_pop=N_pop_true)
    }else{
      dist_dat_temp = data.frame(Sample.Label=i,step_id=k,distance=d_temp_flat,
                                 distance_abs=d_temp_abs,Effort=N_stops,in_vid,
                                 in_angle,Region.Label="Sims",Area=area,N_pop=N_pop_true)
      dist_dat = rbind(dist_dat,dist_dat_temp)
    }
  }
}


dist_dat_in_angle = dist_dat %>% filter(in_angle)



ggplot(dist_dat_in_angle %>% filter(distance<=5),aes(x=distance)) +
  geom_histogram(bins=20) + theme_classic() + xlab("distance (m)")


write.csv(dist_dat_in_angle %>% filter(distance<=5),"results/sim_fig1_dat_alt.csv")



