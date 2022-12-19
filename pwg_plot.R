##version 1.0

####packages####
library(data.table)
library(ggplot2)
library(optparse)
library(jsonlite)

option_list = list(
  make_option(c("-i", "--json"), type="character", default=NULL, help="json file", metavar="character"),
)  

#read json

options(bitmapType='cairo')
svg(paste0(sample_name,".pWGS.svg"), width = 5, height = 1.5)

ggplot(all_results) + 
  geom_jitter(aes(x=0,y=detection_rate,color=type,size=type,shape=type),width = 0.01) +
  
  geom_hline(yintercept = 0,alpha=0.25,color="white") +

  annotate(x = -0.1, xend=0.1, y=dataset_cutoff, yend=dataset_cutoff,
           geom="segment",linetype="dashed",
           colour = "red") +
  
  annotate(geom="text",y = dataset_cutoff,x=0,color="red",label="Detection Cutoff", hjust = 0.5, vjust = -5,size=3) +
  
  #guides(size="none")+
  labs(x="",y="Detection Rate",color="",title="",shape="",size="") +
  scale_color_manual( values= c( "gray", rgb(101/255, 188/255, 69/255) ) ) +
  scale_shape_manual(values=c(1,13)) +
  theme_classic() +
  theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        legend.title=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
 # theme(legend.position="top") +
  coord_flip(clip = "off", xlim=c(-0.1,0.1)) 

dev.off()


