# mrdetect

Workflow for MRdetect, detection of Minimal Residual Disease from paired tumor-plasma sample

## Overview

## Dependencies

* [mrdetect 1.0](https://ctl.cornell.edu/technology/mrdetect-license-request/)
* [bcftools 1.9](https://github.com/samtools/bcftools)


## Usage

### Cromwell
```
java -jar cromwell.jar run mrdetect.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumorSampleName`|String|ID for WGS tumor sample, must match .vcf header
`tumorvcf`|File|tumor vcf file, bgzip
`tumorvcfindex`|File|tumor vcf index file
`reference`|String|genome reference build. Only hg38 supported
`instrument`|String|sequencing instrument used (Illumina NovaSeq X Plus or Illumina NovaSeq 6000)


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`plasmabam`|File?|None|plasma input .bam file
`plasmabai`|File?|None|plasma input .bai file
`plasmaSampleName`|String?|None|name for plasma sample (from bam)
`full_analysis_mode`|Boolean|true|Enable full analysis mode with this flag


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterVCF.tumorVCFfilter`|String|"FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"|set of filter calls to exclude in tumor VCF (any line with these flags will be excluded
`filterVCF.tumorVAF`|String|"0.1"|Variant Allele Frequency for tumor VCF
`filterVCF.jobMemory`|Int|64|Memory allocated for this job (GB)
`filterVCF.threads`|Int|4|Requested CPU threads
`filterVCF.timeout`|Int|10|Hours before task timeout
`parseControls.jobMemory`|Int|4|Memory for this task in GB
`parseControls.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`controlPullReads.jobMemory`|Int|64|Memory allocated for this job (GB)
`controlPullReads.threads`|Int|4|Requested CPU threads
`controlPullReads.timeout`|Int|20|Hours before task timeout
`controlPullReads.pullreadsScript`|String|"$MRDETECT_ROOT/bin/pull_reads"|pull_reads.py executable
`controlQualityScore.jobMemory`|Int|64|Memory allocated for this job (GB)
`controlQualityScore.threads`|Int|4|Requested CPU threads
`controlQualityScore.timeout`|Int|20|Hours before task timeout
`controlQualityScore.pickle`|String|"$MRDETECT_ROOT/bin/MRDetectSNV/trained_SVM.pkl"|trained pickle for detecting real tumor reads
`controlQualityScore.qualityscoreScript`|String|"$MRDETECT_ROOT/bin/quality_score"|quality_score.py executable
`controlDetectSNVs.tumorSampleName`|String|basename(tumorvcf,".vcf")|name for tumour sample (from vcf)
`controlDetectSNVs.jobMemory`|Int|64|Memory allocated for this job (GB)
`controlDetectSNVs.threads`|Int|4|Requested CPU threads
`controlDetectSNVs.timeout`|Int|20|Hours before task timeout
`controlDetectSNVs.blocklist`|String|"$PWGS_BLOCKLIST_ROOT/blocklist.vcf.gz"|list of sites to exclude from analysis, gzipped
`controlDetectSNVs.filterAndDetectScript`|String|"$MRDETECT_ROOT/bin/filterAndDetect"|filterAndDetect.py executable
`samplePullReads.jobMemory`|Int|64|Memory allocated for this job (GB)
`samplePullReads.threads`|Int|4|Requested CPU threads
`samplePullReads.timeout`|Int|20|Hours before task timeout
`samplePullReads.pullreadsScript`|String|"$MRDETECT_ROOT/bin/pull_reads"|pull_reads.py executable
`sampleQualityScore.jobMemory`|Int|64|Memory allocated for this job (GB)
`sampleQualityScore.threads`|Int|4|Requested CPU threads
`sampleQualityScore.timeout`|Int|20|Hours before task timeout
`sampleQualityScore.pickle`|String|"$MRDETECT_ROOT/bin/MRDetectSNV/trained_SVM.pkl"|trained pickle for detecting real tumor reads
`sampleQualityScore.qualityscoreScript`|String|"$MRDETECT_ROOT/bin/quality_score"|quality_score.py executable
`sampleDetectSNVs.tumorSampleName`|String|basename(tumorvcf,".vcf")|name for tumour sample (from vcf)
`sampleDetectSNVs.jobMemory`|Int|64|Memory allocated for this job (GB)
`sampleDetectSNVs.threads`|Int|4|Requested CPU threads
`sampleDetectSNVs.timeout`|Int|20|Hours before task timeout
`sampleDetectSNVs.blocklist`|String|"$PWGS_BLOCKLIST_ROOT/blocklist.vcf.gz"|list of sites to exclude from analysis, gzipped
`sampleDetectSNVs.filterAndDetectScript`|String|"$MRDETECT_ROOT/bin/filterAndDetect"|filterAndDetect.py executable
`snvDetectionSummary.pvalue`|String|"0.00001"|p-value for HBC error rate
`snvDetectionSummary.jobMemory`|Int|20|Memory allocated for this job (GB)
`snvDetectionSummary.threads`|Int|1|Requested CPU threads
`snvDetectionSummary.timeout`|Int|2|Hours before task timeout
`snvDetectionSummary.modules`|String|"mrdetect/1.1.1"|Required environment modules
`snvDetectionSummary.pwgtestscript`|String|"$MRDETECT_ROOT/bin/pwg_test"|executable of pwg_test.R


### Outputs

Output | Type | Description
---|---|---
`snvDetectionResult`|File?|{'description': 'Result from SNV detection incl sample HBCs', 'vidarr_label': 'snvDetectionResult'}
`pWGS_svg`|File?|{'description': 'pWGS svg', 'vidarr_label': 'pWGS_svg'}
`snpcount`|File|{'description': 'number of SNPs in vcf after filtering', 'vidarr_label': 'snpcount'}
`snvDetectionVAF`|File?|{'description': 'VAF from SNV detection for sample', 'vidarr_label': 'snvDetectionVAF'}
`final_call`|File?|{'description': 'Final file of mrdetect results', 'vidarr_label': 'final_call'}
`filteredvcf`|File?|Filtered vcf


## Commands
 This section lists commands run by the MRDetect workflow.
 
 ### filterVCF
 Performs vcf Filtering, followed by processing of individual `MRDetect` calls. Filters include removing difficult regions (optional), splitting multiallelic loci into one allele per line, removing indels, removing loci by quality metrics (set by `tumorVCFfilter`) and finally removing SNPs by VAF (set by `tumorVAF`).
 
 <<<
 		set -euo pipefail
 
 		$BCFTOOLS_ROOT/bin/bcftools view -s ~{tumorSampleName} --regions-file ~{difficultRegions} ~{tumorvcf} |\
 		$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome} |\
 		$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='snps'" |\
 		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{tumorVCFfilter}" |\
 		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{tumorVAF}" > ~{tumorSampleName}.SNP.vcf
 
 		awk '$1 !~ "#" {print}' ~{tumorSampleName}.SNP.vcf | wc -l > ~{tumorSampleName}.SNP.count.txt
 
 	>>>
 
 
 `MRDetect` proceed across three steps. These tasks are run through for the sample and for all the controls.
 
 ### pullReads
 1- `pull_reads` takes any reads in the plasma .bam that corresponds to a SNP in the solid-tumour .vcf. 
 
 <<<
 		set -euo pipefail
 
 		~{pullreadsScript} \
 			--bam ~{plasmabam} \
 			--vcf ~{tumorvcf} \
 			--out PLASMA_VS_TUMOR.tsv
 
 	>>>
 
 ### calculateQualityScore
 2- `quality_score` assesses the likelihood that any read is plasma based on the quality score and the trained pickle. 
 
 <<<
 		set -euo pipefail
 
 		~{qualityscoreScript} \
 			--pickle-name ~{pickle} \
 			--detections ~{snvDetectionReads} \
 			--output_file PLASMA_VS_TUMOR.svm.tsv
 
 	>>>
 
 ### detectSNVs
 3- `filterAndDetect` keeps reads with high plasma likehood and removed those for which SNPs are in the blocklist. 
 
 <<<
 		set -euo pipefail
 
 		~{filterAndDetectScript} \
 			--vcfid ~{tumorSampleName} --bamid ~{plasmaSampleName} \
 			--svm ~{snvDetectionReadsScored} \
 			--vcf ~{tumorvcf} \
 			--output ./ \
 			--blocklist ~{blocklist} \
 	>>>
 
 ### parseControls
 This command processes the list of control files as paired bam/bai files and prints them out for detection.
 
 <<<
 		python <<CODE
 		import os
 
 		with open(os.environ.get("~{controlFileList}")) as f:
 			for line in f:
 				line = line.rstrip()
 				tmp = line.split("\t")
 				r = tmp[0] + "\t" + tmp[1]
 				print(r)
 		f.close()
 		CODE
 	>>>
 
 ### snvDetectionSummary
 
 Finally, `pwg_test.R` will process the controls and the sample to make a final call. 
 
 <<<
 		set -euo pipefail
 
 		cat ~{sep=' ' controlCalls} | awk '$1 !~ "BAM" {print}' > HBCs.csv
 
 		cat ~{sampleCalls} HBCs.csv > ~{plasmaSampleName}.HBCs.csv
 
 		~{pwgtestscript} \
 			--sampleName ~{plasmaSampleName} \
 			--results ~{plasmaSampleName}.HBCs.csv \
 			--candidateSNVsCountFile ~{snpcount} \
 			--vafFile ~{vafFile} \
 			--pval ~{pvalue} 
 
 	>>>
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
