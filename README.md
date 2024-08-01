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


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`plasmabam`|File?|None|plasma input .bam file
`plasmabai`|File?|None|plasma input .bai file
`plasmaSampleName`|String?|None|name for plasma sample (from bam)
`full_analysis_mode`|Boolean|true|Enable full analysis mode with this flag
`controlFileList`|String|"/.mounts/labs/gsi/src/pwgs_hbc/1.0/HBC.bam.list"|tab seperated list of bam and bai files for healthy blood controls


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterVCF.tumorVCFfilter`|String|"FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"|set of filter calls to exclude in tumor VCF (any line with these flags will be excluded
`filterVCF.tumorVAF`|String|"0.1"|Variant Allele Frequency for tumor VCF
`filterVCF.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome .fa
`filterVCF.difficultRegions`|String|"--regions-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"|Path to .bed excluding difficult regions, string must include the flag --regions-file 
`filterVCF.modules`|String|"bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0"|Required environment modules
`filterVCF.jobMemory`|Int|64|Memory allocated for this job (GB)
`filterVCF.threads`|Int|4|Requested CPU threads
`filterVCF.timeout`|Int|10|Hours before task timeout
`parseControls.jobMemory`|Int|4|Memory for this task in GB
`parseControls.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`detectControl.tumorSampleName`|String|basename(tumorvcf,".vcf")|name for tumour sample (from vcf)
`detectControl.modules`|String|"mrdetect/1.1.1 pwgs-blocklist/hg38.1"|Required environment modules
`detectControl.jobMemory`|Int|64|Memory allocated for this job (GB)
`detectControl.threads`|Int|4|Requested CPU threads
`detectControl.timeout`|Int|20|Hours before task timeout
`detectControl.pickle`|String|"$MRDETECT_ROOT/bin/MRDetectSNV/trained_SVM.pkl"|trained pickle for detecting real tumor reads
`detectControl.blocklist`|String|"$PWGS_BLOCKLIST_ROOT/blocklist.vcf.gz"|list of sites to exclude from analysis, gzipped
`detectControl.pullreadsScript`|String|"$MRDETECT_ROOT/bin/pull_reads"|pull_reads.py executable
`detectControl.qualityscoreScript`|String|"$MRDETECT_ROOT/bin/quality_score"|quality_score.py executable
`detectControl.filterAndDetectScript`|String|"$MRDETECT_ROOT/bin/filterAndDetect"|filterAndDetect.py executable
`detectSample.tumorSampleName`|String|basename(tumorvcf,".vcf")|name for tumour sample (from vcf)
`detectSample.modules`|String|"mrdetect/1.1.1 pwgs-blocklist/hg38.1"|Required environment modules
`detectSample.jobMemory`|Int|64|Memory allocated for this job (GB)
`detectSample.threads`|Int|4|Requested CPU threads
`detectSample.timeout`|Int|20|Hours before task timeout
`detectSample.pickle`|String|"$MRDETECT_ROOT/bin/MRDetectSNV/trained_SVM.pkl"|trained pickle for detecting real tumor reads
`detectSample.blocklist`|String|"$PWGS_BLOCKLIST_ROOT/blocklist.vcf.gz"|list of sites to exclude from analysis, gzipped
`detectSample.pullreadsScript`|String|"$MRDETECT_ROOT/bin/pull_reads"|pull_reads.py executable
`detectSample.qualityscoreScript`|String|"$MRDETECT_ROOT/bin/quality_score"|quality_score.py executable
`detectSample.filterAndDetectScript`|String|"$MRDETECT_ROOT/bin/filterAndDetect"|filterAndDetect.py executable
`snvDetectionSummary.pvalue`|String|"0.00001"|p-value for HBC error rate
`snvDetectionSummary.jobMemory`|Int|20|Memory allocated for this job (GB)
`snvDetectionSummary.threads`|Int|1|Requested CPU threads
`snvDetectionSummary.timeout`|Int|2|Hours before task timeout
`snvDetectionSummary.modules`|String|"mrdetect/1.1.1"|Required environment modules
`snvDetectionSummary.pwgtestscript`|String|"$MRDETECT_ROOT/bin/pwg_test"|executable of pwg_test.R


### Outputs

Output | Type | Description | Labels
---|---|---|---
`snvDetectionResult`|File?|Result from SNV detection incl sample HBCs|vidarr_label: snvDetectionResult
`pWGS_svg`|File?|pWGS svg|vidarr_label: pWGS_svg
`snpcount`|File|number of SNPs in vcf after filtering|vidarr_label: snpcount
`snvDetectionVAF`|File?|VAF from SNV detection for sample|vidarr_label: snvDetectionVAF
`final_call`|File?|Final file of mrdetect results|vidarr_label: final_call
`filteredvcf`|File?|Filtered vcf|


## Commands
This section lists commands run by the MRDetect workflow.
 
### filterVCF
Performs vcf Filtering, followed by processing of individual `MRDetect` calls. Filters include removing difficult regions (optional), splitting multiallelic loci into one allele per line, removing indels, removing loci by quality metrics (set by `tumorVCFfilter`) and finally removing SNPs by VAF (set by `tumorVAF`).

```
 
	$BCFTOOLS_ROOT/bin/bcftools view -s ~{tumorbasename} ~{difficultRegions} ~{tumorvcf} |\
	$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome}  |\
	$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='snps'" |\
	$BCFTOOLS_ROOT/bin/bcftools filter -e "~{tumorVCFfilter}" |\
	$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{tumorVAF}" >~{tumorbasename}.SNP.vcf
``` 
### parseControls

This command processes the list of control files as paired bam/bai files and prints them out for detection.
 
### detectSNVs

`MRDetect` proceed across three steps. This task is run through for the sample and for all the controls.
 
* `pull_reads` takes any reads in the plasma .bam that corresponds to a SNP in the solid-tumour .vcf. 
* `quality_score` assesses the likelihood that any read is plasma based on the quality score and the trained pickle. 
* `filterAndDetect` keeps reads with high plasma likehood and removed those for which SNPs are in the blocklist. 
 
```
	$MRDETECT_ROOT/bin/pull_reads \
		--bam ~{plasmabam} \
		--vcf ~{tumorvcf} \
		--out PLASMA_VS_TUMOR.tsv

	$MRDETECT_ROOT/bin/quality_score \
		--pickle-name ~{pickle} \
		--detections PLASMA_VS_TUMOR.tsv \
		--output_file PLASMA_VS_TUMOR.svm.tsv

	$MRDETECT_ROOT/bin/filterAndDetect \
		--vcfid ~{tumorSampleName} --bamid ~{plasmaSampleName} \
		--svm PLASMA_VS_TUMOR.svm.tsv \
		--vcf ~{tumorvcf} \
		--output ./ \
		--blocklist ~{blocklist} \
		--troubleshoot
```
	
### snvDetectionSummary
 
Finally, `pwg_test.R` will process the controls and the sample to make a final call. 	
 		
 		cat ~{sep=' ' controlCalls} | awk '$1 !~ "BAM" {print}' >~{samplebasename}.HBCs.txt
 
 		~{pwgtestscript} \
			--sampleName ~{plasmaSampleName} \
			--results ~{plasmaSampleName}.HBCs.csv \
			--candidateSNVsCountFile ~{snpcount} \
			--vafFile ~{vafFile} \
			--pval ~{pvalue} 

 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
