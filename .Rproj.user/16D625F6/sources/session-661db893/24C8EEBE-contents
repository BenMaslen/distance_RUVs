library(dplyr)

WA_dat = read.csv("data/WA_3DPoints.csv")

#data prep

WA_dat$site_habitat = substr(WA_dat$OpCode,1,3)
WA_dat$site = substr(WA_dat$OpCode,1,2)
WA_dat$habitat = substr(WA_dat$OpCode,3,3)

WA_dat$Family_Genus_Species = as.character(factor(WA_dat$Family):factor(WA_dat$Genus):factor(WA_dat$Species))

WA_dat_fish =  WA_dat[!(WA_dat$Species %in% c("UPPER","LOWER","BOTTOM", "LEFT", "RIGHT")),]

WA_dat_fish_la_ne_po = WA_dat_fish %>% filter(Family_Genus_Species %in% c("Labridae:Choerodon:cyanodus",
                                              "Nemipteridae:Scaevius:milii","Pomacentridae:Dischistodus:darwiniensis"))


#first work out the MeanCount and MaxN per video
pervideo = WA_dat_fish_la_ne_po %>% group_by(Family_Genus_Species,OpCode,site_habitat,FrameRight) %>% 
  count() %>% group_by(Family_Genus_Species,OpCode,site_habitat) %>% 
  summarise(MeanCount=sum(n)/40,MaxN=max(n),sumcount=sum(n))

#then find the number of videos per site
num_vids = WA_dat_fish %>% group_by(site_habitat) %>% 
  summarise(num_vids=length(unique(OpCode)))

pervideo = pervideo %>% left_join(num_vids,by="site_habitat")

#find the average meancount and maxn
maxn_meancount = pervideo %>% group_by(Family_Genus_Species,site_habitat) %>%
  summarise(MeanCount=sum(MeanCount)/max(num_vids),MaxN=sum(MaxN)/max(num_vids),SumCount=sum(sumcount))


write.csv(maxn_meancount,"results/maxn_meancount_WA.csv")

