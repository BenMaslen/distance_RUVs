library(dplyr)
library(Distance)

source("code/RUV_dist_func.R")

WA_dat = read.csv("data/WA_3DPoints.csv")

#data prep
#create unique species and site labels
WA_dat = WA_dat %>% mutate(
  Family_Genus_Species = paste(Family,Genus,Species,sep="_"),
  site_habitat = substr(OpCode,1,3)
)

WA_dat_spatial = WA_dat[WA_dat$Species 
      %in% c("UPPER","LOWER","BOTTOM", "LEFT", "RIGHT"),]
WA_dat_fish =  WA_dat[!(WA_dat$Species 
      %in% c("UPPER","LOWER","BOTTOM", "LEFT", "RIGHT")),]

#Convert into distance data format, estimating 
#distance along the seafloor, horizontal viewing angle
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

#Filter to species Choerodon cyanodus
dist_dat_CC = dist_dat %>% 
  filter(species == "Labridae_Choerodon_cyanodus")

#Left truncation list
trunc.list <- list(left=max(dist_dat_CC$left_trunc,na.rm=T), 
                   right=d_max)

#convert units to hectare
conversion <- convert_units("Metre", NULL, "Hectare")

hr0 <- ds(dist_dat_CC, transect = "point", key="hr", convert_units = conversion, adjustment = NULL,
          formula=~Region.Label,truncation = trunc.list)

density_CC <- dht2(ddf=hr0, flatfile=dist_dat_CC,
                              strat_formula = ~Region.Label, convert_units = conversion,
                              sample_fraction = theta/360,er_est = "P2",stratification = "geographical")



density_CC




