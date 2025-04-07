#Distance sampling simulations

library(tidyverse)
library(MASS)
library(purrr)
library(Distance)
library(doRNG)

set.seed(12345)
overall_population_mean = 100
#pop_disp = 4
N_vids = 10
R_max = 15
theta_diff = 60
gamma_diff = 60
vid_length = 20*60
skipper = 30
N_stops = vid_length/skipper
z_max = 1.5
area = 100*100
N_sims = 1
imp_sd = 5
speed = 10
#file_name = "../../results/Chapter_3/simulation_sims_1000_vids_10_pop_100_2.csv"
schooling=5

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

#set the imperfect detection

# 
norm_x = 0:30

imp_sd = 1


dist_prob_fun <- function(x){
  dnorm(x, mean = 0, sd = imp_sd, log = FALSE)/dnorm(0, mean = 0, sd = imp_sd, log = FALSE)
}


norm_y = dist_prob_fun(norm_x)
imp_func_dat <- data.frame(norm_x,norm_y)

ggplot(imp_func_dat,aes(x=norm_x,y=norm_y)) + geom_line() +
  xlab("Distance") + ylab("Probability of detection") +
  theme_classic()

dist_dat_in_angle$detected = NA

for (l in 1:length(dist_dat_in_angle$distance)){
  dist_dat_in_angle$detected[l] = rbinom(1,1,dist_prob_fun(dist_dat_in_angle$distance[l]))
}

dist_dat_detected <- dist_dat_in_angle  %>% filter(detected==1)

#Histogram of detections by distance
ggplot(dist_dat_detected,aes(x=distance)) +
  geom_histogram(binwidth = 2) + ylab("Detections") +
  xlab("Distance") + theme_classic()


MaxN_dat_prep = dist_dat_detected %>% group_by(step_id,Sample.Label) %>% summarise(MaxNs = sum(detected))
MaxN_dat = MaxN_dat_prep %>% group_by(Sample.Label) %>% summarise(MaxN = max(MaxNs))

MaxN_max  = max(MaxN_dat$MaxN)
MaxN_mean = sum(MaxN_dat$MaxN)/N_vids #do we want the max or the mean of maxes?
MeanCount = sum(dist_dat_detected$detected)/(N_stops*N_vids)


#clean data
dist_dat_detected$Sample.Label <- factor(dist_dat_detected$Sample.Label)
dist_dat_detected <- dist_dat_detected %>% dplyr::select(-detected)

#vids with no fish seen
Empty_vids = (1:N_vids)[!(factor(1:N_vids) %in% levels(dist_dat_detected$Sample.Label))]

#append them into the dataset

#append them into the dataset with catch if all vids have fish
if(length(Empty_vids)!=0){ 
  dist_dat_empty = data.frame(Sample.Label=factor(Empty_vids),step_id=NA,distance=NA,
                              distance_abs=NA,Effort=N_stops,in_vid=FALSE,
                              in_angle=FALSE,Region.Label="Sims",Area=area,N_pop=N_pop_true)
  dist_dat_final = rbind(dist_dat_detected,dist_dat_empty)
}else{
  dist_dat_final = dist_dat_detected
}


#generate object id
dist_dat_final$object <- NA
dist_dat_final$object[!is.na(dist_dat_final$distance)] <- 1:sum(!is.na(dist_dat_final$distance))


#analyse the data

trunc.list <- list(left=z_max/tan(gamma_diff*pi/180), right=R_max)
conversion <- convert_units("meter", NULL, "square meter")

hn0 <- ds(dist_dat_final, transect = "point", key="hn", adjustment = NULL,
          convert_units = conversion, truncation = trunc.list)



# Effective detection radius - need to understand more what this is!

# p_a <- hn0$ddf$fitted[1] #p_a
# w <- R_max-z_max/tan(gamma_diff*pi/180) #range of difference values detected - or the cutpoints
# rho <- sqrt(p_a * w^2) #effective detection radius
# rho


par(mfrow=c(1,2))
plot(hn0, xlab="Distance (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 15)) #detection probability
plot(hn0, xlab="Distance (m)", pdf=TRUE,
     showpoints=FALSE, lwd=3, xlim=c(0, 15)) #Daytime activity


#density estimation
samfrac <- theta_diff / 360


peak.uni.dens <- dht2(hn0, flatfile=dist_dat_final, strat_formula = ~1,
                      sample_fraction = samfrac, er_est = "P2")
# print(peak.uni.dens, report="density")
# 


#create my own bootstrap
nboot = 100

boot_id = rep(NA,nboot)
boot_abund = rep(NA,nboot)
boot_se = rep(NA,nboot)

boot_data = data.frame(boot_id,boot_abund,boot_se)

for (j in 1:nboot){

resample_levels = levels(dist_dat_final$Sample.Label)
boot_sample_levels = sample(resample_levels,replace=T)

boot_dat = data.frame()

for (i in 1:length(boot_sample_levels)){
  samp_level_temp = boot_sample_levels[i]
  boot_dat_level_temp = dist_dat_final %>% filter(Sample.Label==samp_level_temp)
  boot_dat_level_temp$Sample.Label = i
  boot_dat = rbind(boot_dat,boot_dat_level_temp)
}

boot_dat$Sample.Label = factor(boot_dat$Sample.Label)

#generate object id
boot_dat$object <- NA
boot_dat$object[!is.na(boot_dat$distance)] <- 1:sum(!is.na(boot_dat$distance))


hn0_boot <- ds(boot_dat, transect = "point", key="hn", adjustment = NULL,
          convert_units = conversion, truncation = trunc.list)

peak.uni.dens_boot <- dht2(hn0_boot, flatfile=boot_dat, strat_formula = ~1,
                      sample_fraction = samfrac, er_est = "P2")
print(j)
boot_data$boot_id[j] = j
boot_data$boot_se[j] = peak.uni.dens_boot$Abundance_se
boot_data$boot_abund[j] = peak.uni.dens_boot$Abundance
}


boot_data$t = (peak.uni.dens$Abundance - boot_data$boot_abund)/boot_data$boot_se

alpha = 0.05
tUL <- quantile(boot_data$t, probs = c(alpha/2, 1-alpha/2),
                   na.rm=TRUE)

T_CI = c(
  peak.uni.dens$Abundance + tUL[1]*peak.uni.dens$Abundance_se,peak.uni.dens$Abundance + tUL[2]*peak.uni.dens$Abundance_se
)


peak.uni.dens$Abundance + (peak.uni.dens$Abundance -mean(boot_data$boot_abund))







str(peak.uni.dens)

hn0$dht$individuals$D$se
hn0$dht$individuals$D$Estimate
hn0$dht$individuals$D$lcl

peak.uni.dens$
peak.uni.dens$Abundance
peak.uni.dens$LCI
peak.uni.dens$UCI

#get bootstrap confidence intervals
mysummary <- function(ests, fit){
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat = ests$individuals$D$Estimate,
                    Dhat = ests$individuals$D$se))
}


n.cores <- parallel::detectCores()
est.boot <- bootdht(model=hn0, flatfile=dist_dat_final,resample_transect = TRUE,
                    summary_fun=mysummary,sample_fraction = samfrac,
                    convert_units=conversion, nboot=10, cores=n.cores - 1)




peak.uni.dens$Abundance + peak.uni.dens$Abundance - mean(est.boot$D*area)

est.boot

alpha <- 0.05
bootci <- quantile(est.boot$D*area, probs = c(alpha/2, 1-alpha/2),
                    na.rm=TRUE)


par(mfrow=c(1,1))
hist(est.boot$Dhat*area, breaks = 20, 
     xlab="Estimated abundance", main="D-hat estimates bootstraps")
abline(v=quantile(est.boot$Dhat*area, probs = c(0.025,0.975), na.rm=TRUE), lwd=2, lty=3)



peak.uni.dens$Abundance
bootci[1]
bootci[2]
coverage_prob = overall_population_mean<bootci[2] & overall_population_mean>bootci[1]
coverage_prob
bootci_width = bootci[2] - bootci[1]
bootci_width
imp_sd
speed
overall_population_mean
N_vids
num_fish_observed = dim(dist_dat_detected)[1]
num_fish_observed
number_of_frames = N_stops*N_vids
number_of_frames  
