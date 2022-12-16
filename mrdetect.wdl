version 1.0

workflow mrdetect {
	input {
		File plasmabam
		File plasmabai
		String outputFileNamePrefix
		String tumorSampleName
		File tumorvcf
		File tumorvcfindex
		String controlFileList = "/.mounts/labs/gsi/src/pwgs_hbc/1.0/HBC.bam.list"
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		tumorvcf: "tumor vcf file, bgzip"
		tumorvcfindex: "tumor vcf index file"
		outputFileNamePrefix: "Prefix for output file"
		tumorSampleName: "ID for WGS tumor sample, must match .vcf header"
		controlFileList: "tab seperated list of bam and bai files for healthy blood controls"
	}

	call filterVCF {
		input:
		tumorvcf = tumorvcf,
		tumorvcfindex = tumorvcfindex,
		tumorSampleName = tumorSampleName
	}

	call parseControls {
		input:
		controlFileList = controlFileList
	}

	scatter (control in parseControls.controlFiles) {
		call detectSNVs as detectControl {
			input:
			plasmabam = control[0],
			plasmabai = control[1],
			tumorvcf = filterVCF.filteredvcf,
			outputFileNamePrefix = outputFileNamePrefix
		}
	}

	call detectSNVs as detectSample {
		input:
		plasmabam = plasmabam,
		plasmabai = plasmabai,
		tumorvcf = filterVCF.filteredvcf,
		outputFileNamePrefix = outputFileNamePrefix
	}

	call snvDetectionSummary {
		input:
		controlCalls = select_all(detectControl.snvDetectionFinalResult),
		sampleCalls = detectSample.snvDetectionFinalResult,
		outputFileNamePrefix = outputFileNamePrefix 
	}

	meta {
		author: "Felix Beaudry"
		email: "fbeaudry@oicr.on.ca"
		description: "Workflow for MRdetect, detection of Minimal Residual Disease from paired tumor-plasma sample"
		dependencies:
		[
			{
				name: "mrdetect/1.0",
				url: "https://ctl.cornell.edu/technology/mrdetect-license-request/"
			},
			{
				name: "bcftools/1.9",
				url: "https://github.com/samtools/bcftools"
			}
		]
		output_meta: {
			snvDetectionVAF: "VAF from SNV detection for sample",
			snvDetectionHBCResult: "Result from SNV detection incl sample HBCs",
			snpcount: "number of SNPs in vcf after filtering"
		}
	}
	output {
		File? snvDetectionVAF = detectSample.snvDetectionVAF
		File snvDetectionHBCResult = snvDetectionSummary.all_calls
		File snpcount = filterVCF.snpcount
	}
}

task filterVCF {
	input {
		File tumorvcf
		File tumorvcfindex
		String tumorSampleName
		String tumorVCFfilter = "FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"
		String tumorVAF = "0.1"
		String genome = "$HG38_ROOT/hg38_random.fa"
		String difficultRegions = "--regions-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"
		String modules = "bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		tumorvcf: "tumor vcf file, bgzip"
		tumorvcfindex: "tumor vcf index file"
		tumorSampleName: "ID for WGS tumor sample"
		tumorVCFfilter: "set of filter calls to incl. in tumor VCF (any line with these flags will be included"
		tumorVAF: "Variant Allele Frequency for tumor VCF"
		genome: "Path to loaded genome .fa"
		difficultRegions: "Path to .bed excluding difficult regions, string must include the flag --regions-file "
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"

	}

	command <<<
		set -euo pipefail

		$BCFTOOLS_ROOT/bin/bcftools view -s ~{tumorSampleName} ~{difficultRegions} ~{tumorvcf} |\
		$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome} |\
		$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='snps'" |\
		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{tumorVCFfilter}" |\
		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{tumorVAF}" > ~{tumorSampleName}.SNP.vcf

		awk '$1 !~ "#" {print}' ~{tumorSampleName}.SNP.vcf | wc -l >SNP.count.txt

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File filteredvcf = "~{tumorSampleName}.SNP.vcf"
		File snpcount = "SNP.count.txt"
	}

	meta {
		output_meta: {
			filteredvcf: "Filtered vcf",
			snpcount: "Number of SNPs after filtering"
		}
	}
}

task detectSNVs {
	input {
		File plasmabam
		File plasmabai
		String outputFileNamePrefix
		File tumorvcf
		String plasmaSampleName = basename(plasmabam, ".bam")
		String tumorSampleName = basename(tumorvcf, ".vcf")
		String modules = "mrdetect/1.0 pwgs-blocklist/hg38.1"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
		String pickle = "$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/trained_SVM.pkl"
		String blocklist = "$PWGS_BLOCKLIST_ROOT/blocklist.vcf.gz"
		String filterAndDetectScript = "$MRDETECT_ROOT/bin/filterAndDetect"
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		outputFileNamePrefix: "Prefix for output file"
		tumorvcf: "filtered tumor vcf file"
		plasmaSampleName: "name for plasma sample"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		pickle: "trained pickle for detecting real tumor reads"
		blocklist: "list of sites to exclude from analysis, gzipped"
		filterAndDetectScript: "location of filter and detect script"
	}

	command <<<
		set -euo pipefail

		$MRDETECT_ROOT/bin/pull_reads \
			--bam ~{plasmabam} \
			--vcf ~{tumorvcf} \
			--out PLASMA_VS_TUMOR.tsv

		$MRDETECT_ROOT/bin/quality_score \
			--pickle-name ~{pickle} \
			--detections PLASMA_VS_TUMOR.tsv \
			--output_file PLASMA_VS_TUMOR.svm.tsv

		~{filterAndDetectScript} \
			-s PLASMA_VS_TUMOR.svm.tsv \
			-v ~{tumorvcf} \
			-V ~{tumorSampleName} -B ~{plasmaSampleName} \
			-o ./ \
			-b ~{blocklist} \
			-t

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File snvDetectionReadsScored = "PLASMA_VS_TUMOR.svm.tsv" 
		File? snvDetectionFinalResult = "mrdetect.results.csv"
		File snvDetectionVAF = "mrdetect.vaf.txt"
	}

	meta {
		output_meta: {
			snvDetectionReadsScored: "Reads with potential for tumor, with their scores",
			snvDetectionFinalResult: "Final result and call from SNV detection",
			snvDetectionVAF: "Variant Allele Frequencies"
		}
	}
}

task parseControls {
	input {
		String controlFileList 
		Int jobMemory = 4
		Int timeout = 12
	}

	parameter_meta {
		controlFileList: "file with list of control files"
		jobMemory: "Memory for this task in GB"
		timeout: "Timeout in hours, needed to override imposed limits"
	}

	command <<<
		python <<CODE
		import os, re

		with open("~{controlFileList}") as f:
			for line in f:
				line = line.rstrip()
				tmp = line.split("\t")
				r = tmp[0] + "\t" + tmp[1]
				print(r)
		f.close()
		CODE
	>>>

	runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
	}

	output {
		Array[Array[File]] controlFiles = read_tsv(stdout())
	}
}

task snvDetectionSummary {
	input {
		File? sampleCalls
		Array[File] controlCalls
		String outputFileNamePrefix
		String zscoreCutoff = 3.09
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
		String modules = "mrdetect-scripts/1.1"
	}

	parameter_meta {
		sampleCalls: "file of detection rate call for sample"
		controlCalls: "array of file of detection rate calls for HBCs"
		outputFileNamePrefix: "Prefix for output file"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		cat ~{sep=' ' controlCalls} | awk '$1 !~ "BAM" {print}' > HBCs.csv

		cat ~{sampleCalls} HBCs.csv >~{outputFileNamePrefix}.HBCs.csv

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File all_calls = "~{outputFileNamePrefix}.HBCs.csv"
	}

	meta {
		output_meta: {
			all_calls : "HBC mrdetect results"
		}
	}
}
