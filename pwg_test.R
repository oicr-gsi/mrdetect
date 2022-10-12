##version 1.0

####packages####
library(data.table)
library(ggplot2)
library(optparse)
library(jsonlite)

options(scipen=999)

# get options
option_list = list(
  make_option(c("-c", "--controls"), type="character", default=NULL, help="control results file path", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, help="sample results file path", metavar="character"),
  make_option(c("-S", "--sampleName"), type="character", default=NULL, help="sample name", metavar="character"),
  make_option(c("-Z", "--zscoreCutoff"), type="integer", default=1.2, help="Z-score cutoff", metavar="integer")
  
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

sample_path <- opt$sample
control_path <- opt$control
sample_name <- opt$sampleName
zscore_cutoff <- opt$zscoreCutoff

##test
#setwd('/Volumes/cgi/scratch/fbeaudry/plasmaWG/TGL49_0143/')
#sample_path <- 'TGL49_0143_Ct_T_WG_T-92_cfDNA_Input_PLASMA_VS_TUMOR_RESULT.csv'
#control_path <- 'HBCs.txt'
#sample_name <- 'TGL49_0143'

#read files and combine
sample_result <- fread(sample_path,header=FALSE)
control_result <- fread(control_path,header=FALSE)
control_result$V2 <- "CONTROL"
all_results <- rbind.data.frame(sample_result,control_result)[,-7]
names(all_results) <- c("sample","type","sites_checked", "reads_checked", "sites_detected", "detection_rate")

#calculate Z-score for sample
zscore <- (all_results$detection_rate[all_results$type == "TUMOR"] -
              mean(all_results$detection_rate))/
              sd(all_results$detection_rate)

#zscore_cutoff = 1.2 cutoff keeps specificity above 80% empirically in Zviran 2020
if(zscore > zscore_cutoff){
  cancer_detected = "TRUE"
}else if(zscore <= zscore_cutoff){
  cancer_detected = "FALSE"
}else{
  print("Error with z-score")
}

dataset_cutoff <- (zscore_cutoff * sd(all_results$detection_rate)) +  mean(all_results$detection_rate)

detection_multiplier <- all_results$detection_rate[all_results$type == "TUMOR"]/dataset_cutoff

mrdetect_call <- list(zscore,cancer_detected,dataset_cutoff,detection_multiplier)
names(mrdetect_call) <- c("zscore","cancer_detected","dataset_cutoff","detection_multiplier")

HBC_means <- colMeans(control_result[,c("V3","V4","V5","V6"),])
HBC_summary <- list(round(as.numeric(HBC_means[1]),2),round(as.numeric(HBC_means[2]),2),round(as.numeric(HBC_means[3]),2),as.numeric(HBC_means[4]))
names(HBC_summary) <- c("HBC_sites_checked", "HBC_reads_checked", "HBC_sites_detected", "HBC_detection_rate")

sample_results <- sample_result[,c("V3","V4","V5","V6")]
sample_summary <- list(sample_results[[1]],as.numeric(sample_results[[2]]),as.numeric(sample_results[[3]]),as.numeric(sample_results[[4]]))
names(sample_summary) <- c("sites_checked", "reads_checked", "sites_detected", "detection_rate")

all.results <- list(sample_name,sample_summary,HBC_summary,mrdetect_call)
names(all.results) <- c("sample_name","sample_summary","HBC_summary","mrdetect_call")

#convert to JSON and write
ListJSON <- jsonlite::toJSON(all.results,pretty=TRUE,auto_unbox=TRUE)

write(ListJSON,file = paste(sample_name,".mrdetect.json",sep=""))

options(bitmapType='cairo')
svg(paste0(sample_name,".pWGS.svg"), width = 5, height = 1.5)

ggplot(all_results) + 
  geom_jitter(aes(x=0,y=detection_rate,color=type,size=type),width = 0.01) +
  
  geom_hline(yintercept = 0,alpha=0.25,color="white") +

  annotate(x = -0.1, xend=0.25, y=datatset_cutoff, yend=datatset_cutoff,
           geom="segment",linetype="dashed",
           colour = "red") +
  
  annotate(geom="text",y = datatset_cutoff,x=0,color="red",label="Detection Cutoff", hjust = -0.1, vjust = -5,size=3) +
  
  guides(size="none")+
  labs(x="",y="Detection Rate",color="",title="") +
  scale_color_manual( values= c( "gray", rgb(101/255, 188/255, 69/255) ) ) +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
 # theme(legend.position="top") +
  coord_flip(clip = "off", xlim=c(-0.1,0.1)) 

dev.off()


