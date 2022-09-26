library(data.table)
library(tidyverse)
library(ggplot2)
library(optparse)

#
option_list = list(
  make_option(c("-c", "--controls"), type="character", default=NULL, help="control results file path", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, help="sample results file path", metavar="character")
)

# get options
opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

sample_path <- opt$sample
control_path <- opt$control

##test
#setwd(/Volumes/)
#sample_path <- 'cgi/scratch/fbeaudry/plasmaWG/TGL49_0143/TGL49_0143_Ct_T_WG_T-92_cfDNA_Input.filter.deduped.realigned.recalibrated_PLASMA_VS_TUMOR_RESULT.csv'
#control_path <- 'cgi/scratch/fbeaudry/plasmaWG/TGL49_0143/pwg_test.sep26.txt'
#sample_result <- fread(sample_path,header=FALSE)
#control_result <- fread(control_path,header=FALSE)


all_results <- rbind.data.frame(sample_result[,c(V2,V2,V3,V4,V5,V6)],control_result)
names(all_results) <- c("sample","control","sites_checked", "reads_checked", "sites_detected", "detection_rate")
all_results <- all_results %>%  mutate_at(c(3:6), as.numeric)

all_results$type <- "controls"
all_results$type[all_results$sample == "TUMOR"] <- "sample"

Z_high <- (1.3*sd(all_results$detection_rate[all_results$sample != "TUMOR"])) +  mean(all_results$detection_rate[all_results$sample != "TUMOR"])

options(bitmapType='cairo')
svg("pWGS.svg", width = 5, height = 1.5)

ggplot(all_results) + 
  geom_hline(yintercept = 0,alpha=0.25,color="white") +

  annotate(x = -0.1, xend=0.25, y=Z_high, yend=Z_high,
           geom="segment",linetype="dashed",
           colour = "red") +
  
  geom_text(y = Z_high,x=0,color="red",label="Detection Cutoff", hjust = -0.1, vjust = -5,size=3) +
  
  geom_jitter(aes(x=0,y=detection_rate,color=type,size=type),width = 0.01) +
  guides(size="none")+
  theme_bw(base_size=10)+
 # theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))  + 
  labs(x="",y="Detection Rate",color="",title="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
 scale_color_manual(
   values=
     c(
       "gray",
       rgb(101/255, 188/255, 69/255)
       )
   ) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
 # theme(legend.position="top") +
  coord_flip(clip = "off", xlim=c(-0.1,0.1)) 

dev.off()


pzscore <- (all_results$detection_rate[all_results$sample == "TUMOR"] -
mean(all_results$detection_rate[all_results$sample != "TUMOR"]))/
sd(all_results$detection_rate[all_results$sample != "TUMOR"])

if(zscore > 1.2){
  print("Tumor Reads Detected")
}else if(zscore <= 1.2){
  print("No More Tumor Reads Detected Than Expected")
}else{
  print("Error with zscore")
}

