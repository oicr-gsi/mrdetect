##version 1.0

####packages####
library(data.table)
library(optparse)
library(jsonlite)
library(dplyr)

options(scipen=999)

# get options
option_list = list(
  make_option(c("-s", "--sampleName"), type="character", default=NULL, help="sample Name", metavar="character"),
  make_option(c("-r", "--results"), type="character", default=NULL, help="results file path", metavar="character"),
  make_option(c("-c", "--coverage"), type="integer", default=NULL, help="plasma coverage", metavar="character"),
  make_option(c("-C", "--controlCoverageFile"), type="character", default=NULL, help="HBC coverage file", metavar="character"),
  make_option(c("-S", "--candidateSNVsCountFile"), type="character", default=NULL, help="file with the number of SNVs", metavar="character"),
  make_option(c("-p", "--pval"), type="numeric", default=0.001, help="p-value cutoff", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)


sampleName <- opt$sampleName
control_cov_path <- opt$controlCoverageFile
results_path <- opt$results
sample_candidate_SNPs_file <- opt$candidateSNVsCountFile
sample_coverage <- opt$coverage
pval_cutoff <- opt$pval

##test##
#sampleName <- "GLCS_0035_Lu_T_WG.V0.001"
#control_cov_path <- '~/Documents/data/plasmaWG/HBC_coverages.txt'
#results_path <- '~/Documents/data/plasmaWG/lod_seracare/GLCS_0035_Lu_T_WG.V0.001.csv'
#sample_candidate_SNPs_file <- '~/Documents/data/plasmaWG/lod_seracare/GLCS_0035_Lu_T_WG.V0.001.SNP.count.txt'
#sample_coverage <- 30
#pval_cutoff <- 0.001

####read files and combine####
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

tf_estimate <- 1 - ( 1 - (adjusted_detection_rate))^(1/sample_coverage)

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


all.results <- 
  list(
    "sampleName"=sampleName,
     "sample_coverage"=sample_coverage,
     "sample_candidate_SNPs"=sample_candidate_SNPs,
     "tumour_fraction_estimate"=tf_estimate,
     "zscore"=zscore,
     "pvalue"=pvalue,
     "dataset_detection_cutoff"=dataset_cutoff,
     "sites_detected"=results_cov$sites_detected[results_cov$label == "THIS SAMPLE"],
     "mean_noise"=mean(results_cov$noise[results_cov$label != "THIS SAMPLE"]),
     "detection_rate"=results_cov$detection_rate[results_cov$label == "THIS SAMPLE"],
     "false_positive_rate"=mean(results_cov$detection_rate[results_cov$label != "THIS SAMPLE"]),
     "cancer_detected"=cancer_detected
)



#convert to JSON and write
ListJSON <- jsonlite::toJSON(all.results,pretty=TRUE,auto_unbox=TRUE)

write(ListJSON,file = paste(sampleName,".mrdetect.json",sep=""))


