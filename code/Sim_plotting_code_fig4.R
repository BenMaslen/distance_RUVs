library(dplyr)
library(ggplot2)
library(tidyr)


# Set the path to your folder containing distance CSV files
folder_path <- "data/Visibility_sims/"
#folder_path <- "data/batch_test/"

# List all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read and bind all CSV files into one dataframe
combined_df <- csv_files %>%
  lapply(read.csv) %>%
  bind_rows()



Vis_comp = combined_df %>% dplyr::select(MaxN_mean,MeanCount,abundance_estimates,imp_sd)



Vis_comp$imp_sd_fc = Vis_comp$visibility_fc = factor(Vis_comp$imp_sd)
levels(Vis_comp$visibility_fc) = 5:20
Vis_comp$visibility = as.integer(Vis_comp$visibility_fc) + 4

Vis_comp_long = Vis_comp %>% pivot_longer(cols=c("MaxN_mean","MeanCount","abundance_estimates"),
                             names_to = "Method", values_to = "Abundance")

Vis_comp_long$Method = factor(Vis_comp_long$Method)
levels(Vis_comp_long$Method) = c("distance sampling", "MaxN", "MeanCount")

ggplot(Vis_comp_long,
       aes(x=factor(visibility),y=Abundance,colour=visibility)) + geom_boxplot() + 
  theme_classic() + facet_grid(Method~.,scales = "free") +
  scale_y_log10() +xlab("visibility (m)") + ylab("abundance (log scale)") +
  scale_colour_continuous(guide = "none")
  #scale_colour_continuous(name = "Visibility (m)")
  #geom_line(aes(group=1),stat="smooth",method = "gam", formula = y ~ s(x, bs = "cs"),
  #          size = 0.6,colour="blue",
            #linetype ="dashed",
  #          alpha = 0.4)

ggsave("plots/figure_4.png",height=1.1*5.12,width=1.1*5.79*(0.98/0.8))



Vis_comp %>% group_by(visibility) %>% summarise(maxn = mean(MaxN_mean),MeanCount=mean(MeanCount),ds=mean(abundance_estimates))

#MaxN 
#goes from 0.706 to 2.58 as a
2.58/0.706
#3.654391 fold multiplicative increase

#MeanCount 
#goes from 0.0323 to 2.58 as a
0.598/0.0323
#18.51393 fold multiplicative increase



