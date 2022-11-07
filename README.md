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
`plasmabam`|File|plasma input .bam file
`plasmabai`|File|plasma input .bai file
`tumorvcf`|File|tumor vcf file, bgzip
`controlFileList`|String|tab seperated list of bam and bai files for healthy blood controls


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`plasmabasename`|String|basename("~{plasmabam}",".filter.deduped.realigned.recalibrated.bam")|Base name for plasma


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`detectSample.tumorbasename`|String|basename("~{tumorvcf}",".filter.deduped.realigned.recalibrated.mutect2.filtered.vcf.gz")|Base name for tumor
`detectSample.plasmabasename`|String|basename("~{plasmabam}",".filter.deduped.realigned.recalibrated.bam")|Base name for plasma
`detectSample.modules`|String|"mrdetect/1.0 bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0 tabix"|Required environment modules
`detectSample.jobMemory`|Int|64|Memory allocated for this job (GB)
`detectSample.threads`|Int|4|Requested CPU threads
`detectSample.timeout`|Int|10|Hours before task timeout
`detectSample.tumorVCFfilter`|String|"FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"|set of filter calls to incl. in tumor VCF (any line with these flags will be included
`detectSample.tumorVAF`|String|"0.1"|Variant Allele Frequency for tumor VCF
`detectSample.pickle`|String|"$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/trained_SVM.pkl"|trained pickle for detecting real tumor reads
`detectSample.blacklist`|String|"$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/blacklist.txt.gz"|list of sites to exclude from analysis, gzipped
`detectSample.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome .fa
`detectSample.difficultRegions`|String|"--regions-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"|Path to .bed excluding difficult regions, string must include the flag --regions-file 
`detectSample.filterAndDetectScript`|String|"$MRDETECT_ROOT/bin/filterAndDetect"|location of filter and detect script
`parseControls.controlFileListLoc`|String|"~{controlFileList}"|location of file with list of control files, for python intake
`parseControls.jobMemory`|Int|4|Memory for this task in GB
`parseControls.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`detectControl.tumorbasename`|String|basename("~{tumorvcf}",".filter.deduped.realigned.recalibrated.mutect2.filtered.vcf.gz")|Base name for tumor
`detectControl.plasmabasename`|String|basename("~{plasmabam}",".filter.deduped.realigned.recalibrated.bam")|Base name for plasma
`detectControl.modules`|String|"mrdetect/1.0 bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0 tabix"|Required environment modules
`detectControl.jobMemory`|Int|64|Memory allocated for this job (GB)
`detectControl.threads`|Int|4|Requested CPU threads
`detectControl.timeout`|Int|10|Hours before task timeout
`detectControl.tumorVCFfilter`|String|"FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"|set of filter calls to incl. in tumor VCF (any line with these flags will be included
`detectControl.tumorVAF`|String|"0.1"|Variant Allele Frequency for tumor VCF
`detectControl.pickle`|String|"$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/trained_SVM.pkl"|trained pickle for detecting real tumor reads
`detectControl.blacklist`|String|"$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/blacklist.txt.gz"|list of sites to exclude from analysis, gzipped
`detectControl.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome .fa
`detectControl.difficultRegions`|String|"--regions-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"|Path to .bed excluding difficult regions, string must include the flag --regions-file 
`detectControl.filterAndDetectScript`|String|"$MRDETECT_ROOT/bin/filterAndDetect"|location of filter and detect script
`snvDetectionSummary.DetectionRScript`|String|"$MRDETECT_SCRIPTS_ROOT/bin/pwg_test.R"|location of pwg_test.R
`snvDetectionSummary.samplebasename`|String|basename("~{sampleCalls}",".PLASMA_VS_TUMOR_RESULT.csv")|base name for files
`snvDetectionSummary.jobMemory`|Int|20|Memory allocated for this job (GB)
`snvDetectionSummary.threads`|Int|1|Requested CPU threads
`snvDetectionSummary.timeout`|Int|2|Hours before task timeout
`snvDetectionSummary.modules`|String|"mrdetect-scripts/1.0"|Required environment modules


### Outputs

Output | Type | Description
---|---|---
`snvDetectionFinalResult`|File|Final result and call from SNV detection
`snvDetectionHBCResult`|File|results from the HBCs
`pWGS_svg`|File|pWGS svg


## Commands
 This section lists command(s) run by WORKFLOW workflow
 
 * Running WORKFLOW
 
 === Description here ===.
 
 <<<
 		set -euo pipefail
 
 		tabix -fp vcf ~{tumorvcf}
 
 		$BCFTOOLS_ROOT/bin/bcftools view -s ~{tumorbasename} ~{difficultRegions} ~{tumorvcf} |\
 		$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome}  |\
 		$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='snps'" |\
 		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{tumorVCFfilter}" |\
 		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{tumorVAF}" >~{tumorbasename}.SNP.vcf
 
 		$MRDETECT_ROOT/bin/pull_reads \
 			--bam ~{plasmabam} \
 			--vcf ~{tumorbasename}.SNP.vcf \
 			--out ~{plasmabasename}_PLASMA_VS_TUMOR
 
 		$MRDETECT_ROOT/bin/quality_score \
 			--pickle-name ~{pickle} \
 			--detections ~{plasmabasename}_PLASMA_VS_TUMOR.tsv \
 			--output_file ~{plasmabasename}_PLASMA_VS_TUMOR.svm.tsv
 
 		cp ~{blacklist} ./blacklist.txt.gz
 
 		~{filterAndDetectScript} \
 			~{tumorbasename}.SNP.vcf \
 			~{plasmabasename}_PLASMA_VS_TUMOR.svm.tsv \
 			~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv >detection_output.txt
 
 		awk '$1 ~ "chr" {print $1"\t"$2"\t"$3"\t"$4}' detection_output.txt | uniq -c >detectionsPerSite.txt
 
 
 	>>>
 <<<
 		python <<CODE
 		import os, re
 
 		with open("~{controlFileListLoc}") as f:
 			for line in f:
 				line = line.rstrip()
 				tmp = line.split("\t")
 				r = tmp[0] + "\t" + tmp[1]
 				print(r)
 		f.close()
 		CODE
 	>>>
 <<<
 		set -euo pipefail
 
 		cat ~{sep=' ' controlCalls} | awk '$1 !~ "BAM" {print}' >~{samplebasename}.HBCs.txt
 
 		Rscript --vanilla ~{DetectionRScript} -s ~{sampleCalls} -S ~{samplebasename} -c ~{samplebasename}.HBCs.txt
 
 	>>>
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
