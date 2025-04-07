library(dplyr)
library(Distance)

source("code/RUV_dist_func.R")

WA_dat = read.csv("data/WA_3DPoints.csv")

#data prep
WA_dat$site_habitat = substr(WA_dat$OpCode,1,3)
WA_dat$Family_Genus_Species = paste(WA_dat$Family,WA_dat$Genus,WA_dat$Species,sep="_")


WA_dat_spatial = WA_dat[WA_dat$Species 
                        %in% c("UPPER","LOWER","BOTTOM", "LEFT", "RIGHT"),]
WA_dat_fish =  WA_dat[!(WA_dat$Species 
                        %in% c("UPPER","LOWER","BOTTOM", "LEFT", "RIGHT")),]


dist_dat = RUVdist.prep(
  spatial_dat = WA_dat_spatial,
  fish_dat = WA_dat_fish,
  vid_name = "OpCode",
  species = "Family_Genus_Species",
  site = "site_habitat",
  cam_height = 350,
  height_scaler = 1.1,
  nframes = 40)

#maximum viewing distance
d_max = 5

#get theta as the smallest of the estimated thetas
theta = min(dist_dat$theta,na.rm = T)

#unit conversion
conversion <- convert_units("Metre", NULL, "Hectare")


########
#filter to species Labridae_Choerodon_cyanodus
WA_dat_distance_labridae = dist_dat %>% filter(species == "Labridae_Choerodon_cyanodus")

#look at distances
hist(WA_dat_distance_labridae$distance)
max(WA_dat_distance_labridae$distance,na.rm=T)

#left and right truncations
trunc.list <- list(left=max(WA_dat_distance_labridae$left_trunc,na.rm=T), right=d_max)


#fit models
hn0 <- ds(WA_dat_distance_labridae, transect = "point", key="hn", convert_units = conversion, adjustment = NULL,
          formula=~Region.Label,truncation = trunc.list)
hn1 <- ds(WA_dat_distance_labridae, transect = "point", key="hn", convert_units = conversion, adjustment = "cos",
          nadj=1,formula=~Region.Label,truncation = trunc.list)
hn2 <- ds(WA_dat_distance_labridae, transect = "point", key="hn", convert_units = conversion, adjustment = "cos",
          nadj=2,formula=~Region.Label,truncation = trunc.list)


QAIC(hn0,hn1)
#choose hn0

hr0 <- ds(WA_dat_distance_labridae, transect = "point", key="hr", convert_units = conversion, adjustment = NULL,
          formula=~Region.Label,truncation = trunc.list)
hr1 <- ds(WA_dat_distance_labridae, transect = "point", key="hr", convert_units = conversion, adjustment = "poly",
          nadj=1,formula=~Region.Label,truncation = trunc.list)

QAIC(hr0,hr1)
#choose hr0

#now use c values to choose between the models
chats <- chi2_select(hn0, hr0)$criteria
modnames <- unlist(lapply(list(hn0, hr0), function(x) x$ddf$name.message))
results <- data.frame(modnames, chats)
results.sort <- results[order(results$chats),]
results.sort

#choose hazard ratio with no adjustment

chosen_model_labridae = hr0

#check
ddf.gof(chosen_model_labridae$ddf)


#look at help(QAIC) for paper to cite on model selection method = How et al (2019), also howe et al 2017 from the vingete)

#model summary
summary(chosen_model_labridae)

#density estimates
WA_fish.ests_labridae <- dht2(ddf=chosen_model_labridae, flatfile=WA_dat_distance_labridae,
                              strat_formula = ~Region.Label, convert_units = conversion,
                              sample_fraction = theta/360,er_est = "P2",stratification = "geographical")



WA_fish.ests_labridae

#bootstrap confidence interval
mysummary <- function(ests, fit){
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat = ests$individuals$D$Estimate))
}

set.seed(12345)

n.cores <- parallel::detectCores()
daytime.boot.uni <- bootdht(model=chosen_model_labridae, flatfile=WA_dat_distance_labridae,
                            resample_transects = TRUE, nboot=1000, 
                            cores = n.cores - 1,
                            summary_fun=mysummary, sample_fraction = theta/360,
                            convert_units = conversion)

print(summary(daytime.boot.uni))

WA_fish.boot_labridae <-  daytime.boot.uni %>% group_by(Label) %>% summarise(LCL=quantile(Dhat, probs = c(0.025), na.rm=TRUE),
                                                                             UCL=quantile(Dhat, probs = c(0.975), na.rm=TRUE))

#write.csv(WA_fish.boot_labridae,"results/bootstrap_lab.csv")

#now lets repeat for the other two species

########
#filter to species Nemipteridae_Scaevius_milii

WA_dat_distance_milli = dist_dat %>% filter(species == "Nemipteridae_Scaevius_milii")

#look at distances
hist(WA_dat_distance_milli$distance)
max(WA_dat_distance_milli$distance,na.rm=T)

#left and right truncations
trunc.list <- list(left=max(WA_dat_distance_milli$left_trunc,na.rm=T), right=d_max)

#remove SSS
WA_dat_distance_milli %>% group_by(Region.Label) %>% summarise(count=sum(distance>0,na.rm=T))

#remove regions with no observations
WA_dat_distance_milli_SSS <- WA_dat_distance_milli

#remove regions with no & small observations (<5)
WA_dat_distance_milli <- WA_dat_distance_milli %>% filter(Region.Label!="SSS") 



#start model selection on model with small 
#observation sites removed as we do not have enough information to estimate the detection function here


hn0 <- ds(WA_dat_distance_milli, transect = "point", key="hn", convert_units = conversion, adjustment = NULL,
          formula=~Region.Label,truncation = trunc.list)
hn1 <- ds(WA_dat_distance_milli, transect = "point", key="hn", convert_units = conversion, adjustment = "cos",
          nadj=1,formula=~Region.Label,truncation = trunc.list)
hn2 <- ds(WA_dat_distance_milli, transect = "point", key="hn", convert_units = conversion, adjustment = "cos",
          nadj=2,formula=~Region.Label,truncation = trunc.list)
QAIC(hn0,hn1,hn2)
#choose hn0 even though it has higher QAIC as it converged properly unlike hn1 and hn2

hr0 <- ds(WA_dat_distance_milli, transect = "point", key="hr", convert_units = conversion, adjustment = NULL,
          formula=~Region.Label,truncation = trunc.list)
hr1 <- ds(WA_dat_distance_milli, transect = "point", key="hr", convert_units = conversion, adjustment = "poly",
          nadj=1,formula=~Region.Label,truncation = trunc.list)

QAIC(hr0,hr1)
#choose hr0



#now use c values to choose between the models
chats <- chi2_select( hn0, hr0)$criteria
modnames <- unlist(lapply(list( hn0, hr0), function(x) x$ddf$name.message))
results <- data.frame(modnames, chats)
results.sort <- results[order(results$chats),]
results.sort

#choose hazard ratio with no adjustment

chosen_model_milli = hr0

#check
ddf.gof(chosen_model_milli$ddf)

#model summary
summary(chosen_model_milli)

#density estimates
WA_fish.ests_milli <- dht2(ddf=chosen_model_milli, flatfile=WA_dat_distance_milli,
                           strat_formula = ~Region.Label, convert_units = conversion,
                           sample_fraction = theta/360,er_est = "P2",stratification = "geographical")


WA_fish.ests_milli

#bootstrap confidence interval
mysummary <- function(ests, fit){
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat = ests$individuals$D$Estimate))
}

n.cores <- parallel::detectCores()
daytime.boot.uni.milli <- bootdht(model=chosen_model_milli, flatfile=WA_dat_distance_milli,
                                  resample_transects = TRUE, nboot=1000, 
                                  cores = n.cores - 1,
                                  summary_fun=mysummary, sample_fraction = theta/360,
                                  convert_units = conversion)

print(summary(daytime.boot.uni.milli))

WA_fish.boot_milli <- daytime.boot.uni.milli %>% group_by(Label) %>% summarise(LCL=quantile(Dhat, probs = c(0.025), na.rm=TRUE),
                                                                               UCL=quantile(Dhat, probs = c(0.975), na.rm=TRUE))



#Estimate SSS

chosen_model_milli_SSS <- ds(WA_dat_distance_milli_SSS, transect = "point", key="hr", convert_units = conversion, adjustment = NULL,
                             truncation = trunc.list)


WA_fish.ests_milli_SSS <- dht2(ddf=chosen_model_milli_SSS, flatfile=WA_dat_distance_milli_SSS,
                               strat_formula = ~Region.Label, convert_units = conversion,
                               sample_fraction = theta/360,er_est = "P2",stratification = "geographical")


WA_fish.ests_milli_SSS <- data.frame(WA_fish.ests_milli_SSS) %>% filter(Region.Label=="SSS")


# bootstrap
daytime.boot.uni.milli_SSS <- bootdht(model=chosen_model_milli_SSS, flatfile=WA_dat_distance_milli_SSS,
                                      resample_transects = TRUE, nboot=1000, 
                                      cores = n.cores - 1,
                                      summary_fun=mysummary, sample_fraction = 51.99647/360,
                                      convert_units = conversion)

print(summary(daytime.boot.uni.milli))

WA_fish.boot_milli_SSS <- daytime.boot.uni.milli_SSS %>% group_by(Label) %>% 
  summarise(LCL=quantile(Dhat, probs = c(0.025), na.rm=TRUE),
            UCL=quantile(Dhat, probs = c(0.975), na.rm=TRUE)) %>%
  filter(Label=="SSS")


########
#filter to species Pomacentridae_Dischistodus_darwiniensis
WA_dat_distance_disch = dist_dat %>% filter(species == "Pomacentridae_Dischistodus_darwiniensis")

#look at distances
hist(WA_dat_distance_disch$distance)
max(WA_dat_distance_disch$distance,na.rm=T)

#left and right truncations
trunc.list <- list(left=max(WA_dat_distance_disch$left_trunc,na.rm=T), right=d_max)

WA_dat_distance_disch %>% group_by(Region.Label) %>% summarise(count=sum(distance>0,na.rm=T))


#start model selection

hn0 <- ds(WA_dat_distance_disch, transect = "point", key="hn", convert_units = conversion, adjustment = NULL,
          formula=~Region.Label,truncation = trunc.list)
hn1 <- ds(WA_dat_distance_disch, transect = "point", key="hn", convert_units = conversion, adjustment = "cos",
          nadj=1,formula=~Region.Label,truncation = trunc.list)
hn2 <- ds(WA_dat_distance_disch, transect = "point", key="hn", convert_units = conversion, adjustment = "cos",
          nadj=2,formula=~Region.Label,truncation = trunc.list)
QAIC(hn0,hn1,hn2)
#choose hn0 (as it converged without issue)

hr0 <- ds(WA_dat_distance_disch, transect = "point", key="hr", convert_units = conversion, adjustment = NULL,
          formula=~Region.Label,truncation = trunc.list)
hr1 <- ds(WA_dat_distance_disch, transect = "point", key="hr", convert_units = conversion, adjustment = "poly",
          nadj=1,formula=~Region.Label,truncation = trunc.list)

QAIC(hr0,hr1)
#choose hr0

#now use c values to choose between the models
chats <- chi2_select(hn0, hr0)$criteria
modnames <- unlist(lapply(list(hn0, hr0), function(x) x$ddf$name.message))
results <- data.frame(modnames, chats)
results.sort <- results[order(results$chats),]
results.sort

#choose hazard rate with no adjustment

chosen_model_disch = hr0

#check
ddf.gof(chosen_model_disch$ddf)

#model summary
summary(chosen_model_disch)

#density estimates
WA_fish.ests_disch <- dht2(ddf=chosen_model_disch, flatfile=WA_dat_distance_disch,
                           strat_formula = ~Region.Label, convert_units = conversion,
                           sample_fraction = theta/360,er_est = "P2",stratification = "geographical")


WA_fish.ests_disch

#bootstrap confidence interval
mysummary <- function(ests, fit){
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat = ests$individuals$D$Estimate))
}

n.cores <- parallel::detectCores()
daytime.boot.uni.disch <- bootdht(model=chosen_model_disch, flatfile=WA_dat_distance_disch,
                                  resample_transects = TRUE, nboot=1000, 
                                  cores = n.cores - 1,
                                  summary_fun=mysummary, sample_fraction = theta/360,
                                  convert_units = conversion)

print(summary(daytime.boot.uni.disch))

WA_fish.boot_disch <- daytime.boot.uni.disch %>% group_by(Label) %>% summarise(LCL=quantile(Dhat, probs = c(0.025), na.rm=TRUE),
                                                                               UCL=quantile(Dhat, probs = c(0.975), na.rm=TRUE))




#join the results together

boot_results = rbind(WA_fish.boot_labridae,WA_fish.boot_milli_SSS,WA_fish.boot_milli,WA_fish.boot_disch)


results_la_ma_di <- rbind(data.frame(WA_fish.ests_labridae),WA_fish.ests_milli_SSS,data.frame(WA_fish.ests_milli),data.frame(WA_fish.ests_disch))
results_la_ma_di$species = rep(c("Labridae_Choerodon_cyanodus",
                                 "Nemipteridae_Scaevius_milii","Pomacentridae_Dischistodus_darwiniensis"),c(dim(data.frame(WA_fish.ests_labridae))[1],
                                                                                                            dim(data.frame(WA_fish.ests_milli))[1]+1,
                                                                                                            dim(data.frame(WA_fish.ests_disch))[1]))

results_la_ma_di <- cbind(results_la_ma_di,boot_results)

#effective detection radius (as obtained from https://examples.distancesampling.org/Distance-cameratraps/camera-distill.html)
p_a_l <- unique(chosen_model_labridae$ddf$fitted)
p_a_m <- c(unique(chosen_model_milli_SSS$ddf$fitted),unique(chosen_model_milli$ddf$fitted))
p_a_d <- unique(chosen_model_disch$ddf$fitted)

w_l <- 5-max(WA_dat_distance_labridae$left_trunc,na.rm=T)
w_m <- 5-max(WA_dat_distance_milli$left_trunc,na.rm=T)
w_d <- 5-max(WA_dat_distance_disch$left_trunc,na.rm=T)

rho_l <- sqrt(p_a_l * w_l^2)
rho_m <- sqrt(p_a_m * w_m^2)
rho_p <- sqrt(p_a_d * w_d^2)

results_la_ma_di$effective_detect_radius = c(rho_l,NA,rho_m,NA,rho_p,NA)



write.csv(results_la_ma_di,"results/results_la_ma_di.csv")

#combine the plots
dev.off()
par(mfrow=c(1,3))


#png(filename="plots/WA_dist_estimates.png")

plot(chosen_model_labridae, showpoints=F, main="Choerodon cyanodus",xlab="",ylab="detection probability",
     lty=0,pl.col = "darkgrey",ylim=c(0,1.25))
add.df.covar.line(chosen_model_labridae, data=data.frame(Region.Label="JYA"), 
                  lwd=3, lty=1, col="#1B9E77")
add.df.covar.line(chosen_model_labridae, data=data.frame(Region.Label="MOA"), 
                  lwd=3, lty=1, col="#7570B3")
add.df.covar.line(chosen_model_labridae, data=data.frame(Region.Label="MOS"), 
                  lwd=3, lty=1, col="#E7298A")
add.df.covar.line(chosen_model_labridae, data=data.frame(Region.Label="SSC"), 
                  lwd=3, lty=1, col="#D95F02")


plot(chosen_model_milli_SSS, showpoints=FALSE, main="Scaevius milii",
     ylab="",xlab="distance (m)",lty=0,pl.col = "darkgrey",ylim=c(0,1.25))
add.df.covar.line(chosen_model_milli, data=data.frame(Region.Label="MOA"), 
                  lwd=3, lty=1, col="#7570B3")
add.df.covar.line(chosen_model_milli, data=data.frame(Region.Label="MOS"), 
                  lwd=3, lty=1, col="#E7298A")
add.df.covar.line(chosen_model_milli_SSS, data=data.frame(Region.Label="SSS"), 
                  lwd=3, lty=1, col="#E6AB02")


plot(chosen_model_disch, showpoints=FALSE, main="Dischistodus species",
     ylab="",xlab="",lty=0,pl.col = "darkgrey",ylim=c(0,1.25))
add.df.covar.line(chosen_model_disch, data=data.frame(Region.Label="JYA"), 
                  lwd=3, lty=1, col="#1B9E77")
add.df.covar.line(chosen_model_disch, data=data.frame(Region.Label="MOA"), 
                  lwd=3, lty=1, col="#7570B3")
add.df.covar.line(chosen_model_disch, data=data.frame(Region.Label="SSC"), 
                  lwd=3, lty=1, col="#D95F02")





par(mfrow=c(1,1))

# Adjust the margins
par(mar=c(5, 4, 4, 2) + 4,xpd=T)


plot(chosen_model_disch, showpoints=FALSE, main="Pomacentridae Dischistodus species",ylab="",xlab="")
legend("bottom", legend=c("Macroalgae - JY", "Macroalgae - MO", "Seagrass - MO","Seagrass - SSI","Coral - SSI"),horiz=T,
       lwd=3, lty=1, inset=c(0,-0.5),col=c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02","#D95F02"))





