library(data.table)
library(tidyverse)
library(cowplot)
library(ggdist)
library(gghalves)

#https://z3tt.github.io/Rainclouds/

setwd('/Volumes/')

MRDetect_raw <- fread('cgi/scratch/fbeaudry/plasmaWG/pwg_test.may25.txt',header=TRUE)

MRDetect_raw$type <- "mismatch"
MRDetect_raw$type[MRDetect_raw$tumor == MRDetect_raw$plasma] <- "match"
 
MRDetect_raw$sitesCheckedMatch <- NA
MRDetect_raw$sitesCheckedMatch[MRDetect_raw$tumor == MRDetect_raw$plasma] <- MRDetect_raw$sitesChecked[MRDetect_raw$tumor == MRDetect_raw$plasma]


detection_rate_plot <- 
ggplot(MRDetect_raw) + 
  geom_hline(yintercept = 0,alpha=0.25) +
  geom_text(aes(x=tumor,y=-0.001,label=sitesCheckedMatch)) +
  
  
#  geom_boxplot(aes(x=tumor,y=detectionRate,color=type),width=0.1) +
  geom_jitter(aes(x=tumor,y=detectionRate,color=type),width = 0.1,alpha=0.75,size=2) +
  
 # ggdist::stat_pointinterval(aes(x=tumor,y=detectionRate,color=type)) +
  
  theme_bw(base_size=15)+
  theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0))  + 
  labs(x="Samples",y="Detection Rate",color="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
 scale_color_manual(
   values=
     c(
       rgb(101/255, 188/255, 69/255), 
       rgb(0,0,0)
       )
   ) +
  guides(color="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

 
MRDetect <- MRDetect_raw %>% group_by(tumor) %>% 
  mutate(zscore=(detectionRate-mean(detectionRate))/sd(detectionRate))

zscore_plot <- 
ggplot(MRDetect %>% filter(type == "match"),aes(x=tumor,y=zscore)) + 
  geom_point(size=3,shape=1) + theme_bw(base_size=15) +   
  geom_hline(yintercept = 1.2,alpha=0.25,linetype="dashed")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0))  +
  scale_color_manual(
    values=
      c(
        rgb(101/255, 188/255, 69/255), 
        rgb(0,0,0)
      )
  ) +   labs(x="Samples",y="Z-score",color="") 

plot_grid(detection_rate_plot,zscore_plot,  ncol = 1, rel_heights = c(1,0.75),align = 'v')


