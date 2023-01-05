##version 1.0

####packages####
library(data.table)
library(optparse)
library(dplyr)

options(scipen=999)
options(digits = 5)

# get options
option_list = list(
  make_option(c("-s", "--sampleName"), type="character", default=NULL, help="sample Name", metavar="character"),
  make_option(c("-r", "--results"), type="character", default=NULL, help="results file path", metavar="character"),
  make_option(c("-v", "--vafFile"), type="character", default=NULL, help="VAF file", metavar="character"),
  make_option(c("-C", "--controlCoverageFile"), type="character", default=NULL, help="HBC coverage file", metavar="character"),
  make_option(c("-S", "--candidateSNVsCountFile"), type="character", default=NULL, help="file with the number of SNVs", metavar="character"),
  make_option(c("-p", "--pval"), type="numeric", default=0.001, help="p-value cutoff", metavar="numeric"),
  make_option(c("-j", "--json"), type="character", default=FALSE, help="export result as json", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

sampleName <- opt$sampleName
control_cov_path <- opt$controlCoverageFile
results_path <- opt$results
sample_candidate_SNPs_file <- opt$candidateSNVsCountFile
pval_cutoff <- opt$pval
as.json <- opt$json
vaf_file <- opt$vafFile

##test##
#sampleName <- "GLCS_0035_Lu_T_WG"
#control_cov_path <- '/Volumes/cgi/scratch/fbeaudry/plasmaWG/HBC_coverages.txt'
#results_path <- '/Volumes/cgi/scratch/fbeaudry/plasmaWG/seracare_lod/insilico/GLCS_0035_Lu_T_WG.HBCs.csv'
#sample_candidate_SNPs_file <- '/Volumes/cgi/scratch/fbeaudry/plasmaWG/seracare_lod/insilico/GLCS_0035_Lu_T_WG.SNP.count.txt'
#pval_cutoff <- 0.001
#as.json <- FALSE
#vaf_file <- '/Volumes/cgi/scratch/fbeaudry/plasmaWG/seracare_lod/insilico/GLCS_0035_Lu_T_WG.mrdetect.vaf.txt'

####read files and combine####
vaf <- fread(vaf_file)
sample_coverage = median(vaf$goodreads)
median_vaf = median(vaf$vaf)

sample_candidate_SNPs <- as.numeric(unlist(fread(sample_candidate_SNPs_file))[[1]])

results <- fread(results_path,header=TRUE)
results$label <- "CONTROLS"
results$label[1] <- "THIS SAMPLE"

hbc_coverage <- fread(control_cov_path)
results_cov <- left_join(results,hbc_coverage)

results_cov$coverage[1] <- sample_coverage

####calculations####
#noise rate = sites detected * coverage
results_cov$noise_rate <- results_cov$sites_detected / results_cov$coverage
results_cov$noise <- results_cov$noise_rate * sample_coverage

#calculate Tumour Fraction (TF) for sample
adjusted_detection_rate <- (results_cov$sites_detected[results_cov$label == "THIS SAMPLE"] -
                              mean(results_cov$noise[results_cov$label != "THIS SAMPLE"])) / 
                              sample_candidate_SNPs

tf_estimate <- 
  1 - ( 1 - (adjusted_detection_rate))^(1/sample_coverage)

#calculate Z-score for sample
zscore <- (results_cov$noise[results_cov$label == "THIS SAMPLE"] -
             mean(results_cov$noise))/
  sd(results_cov$noise)

pvalue <- pnorm(zscore,lower.tail=F)


if(pvalue < pval_cutoff){
  cancer_detected = "TRUE"
}else if(pvalue >= pval_cutoff){
  cancer_detected = "FALSE"
}else{
  print("Error with p-value")
}

dataset_cutoff <- (qnorm(pval_cutoff,lower.tail = F) * sd(results_cov$detection_rate)) +  mean(results_cov$detection_rate)

if(as.json == TRUE){
  library(jsonlite)
  
  all.results <- 
    list(
      "sampleName"=sampleName,
      "sample_coverage"=sample_coverage,
      "median_vaf"=median_vaf,
      "sample_candidate_SNPs"=sample_candidate_SNPs,
      "sites_detected"=results_cov$sites_detected[results_cov$label == "THIS SAMPLE"],
      "mean_noise"=mean(results_cov$noise[results_cov$label != "THIS SAMPLE"]),
      "detection_rate"=results_cov$detection_rate[results_cov$label == "THIS SAMPLE"],
      "tumour_fraction_estimate"=tf_estimate,
      "zscore"=zscore,
      "pvalue"=pvalue,
      "dataset_detection_cutoff"=dataset_cutoff,
      "false_positive_rate"=mean(results_cov$detection_rate[results_cov$label != "THIS SAMPLE"]),
      "cancer_detected"=cancer_detected
      )
  
  #convert to JSON and write
  ListJSON <- jsonlite::toJSON(all.results,pretty=TRUE,auto_unbox=TRUE)
  
  write(ListJSON,file = paste(sampleName,".mrdetect.json",sep=""))
} else {
 results.table <- cbind(
    "sampleName"=sampleName,
    "sample_coverage"=sample_coverage,
    "median_vaf"=median_vaf,
    "sample_candidate_SNPs"=sample_candidate_SNPs,
    "sites_detected"=results_cov$sites_detected[results_cov$label == "THIS SAMPLE"],
    "mean_noise"=mean(results_cov$noise[results_cov$label != "THIS SAMPLE"]),
    "detection_rate"=results_cov$detection_rate[results_cov$label == "THIS SAMPLE"],
    "tumour_fraction_estimate"=tf_estimate,
    "zscore"=zscore,
    "pvalue"=pvalue,
    "dataset_detection_cutoff"=dataset_cutoff,
    "false_positive_rate"=mean(results_cov$detection_rate[results_cov$label != "THIS SAMPLE"]),
    "cancer_detected"=cancer_detected
  )
 write.table(
   results.table,
   file = paste(sampleName,".mrdetect.txt",sep=""),
   append = F, quote = FALSE, sep = "\t", 
   eol = "\n", na = "NA",dec = ".", row.names = FALSE, 
   col.names = TRUE
 )
}

