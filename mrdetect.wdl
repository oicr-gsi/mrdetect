version 1.0

workflow mrdetect {
	input {
		File plasmabam 
		File plasmabai
		String plasmabasename = basename("~{plasmabam}", ".filter.deduped.realigned.recalibrated.bam")
		File tumorvcf
		File controlFileList
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		tumorvcf: "tumor vcf file, bgzip"
		plasmabasename: "Base name for plasma"
		controlFileList: "tab seperated list of bam and bai files for healthy blood controls"
	}

	call detectSNVs as detectSample {
		input: 
		plasmabam = plasmabam, 
		plasmabai = plasmabai, 
		tumorvcf = tumorvcf
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
			tumorvcf = tumorvcf
		}
	}

	call snvDetectionSummary {
		input:
		controlCalls = select_all(detectControl.snvDetectionFinalResult),
		sampleCalls = detectSample.snvDetectionFinalResult
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
			snvDetectionFinalResult: "Final result and call from SNV detection",
			pWGS_svg: "pWGS svg",
			snvDetectionHBCResult: "results from the HBCs"
		}
	}
	output {
		File snvDetectionFinalResult = "~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv"
		File snvDetectionHBCResult = "~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv.HBCs.txt"
		File pWGS_svg = "pWGS.svg"
	}
}

task detectSNVs {
	input {
		File plasmabam 
		File plasmabai 
		File tumorvcf
		String tumorbasename = basename("~{tumorvcf}", ".filter.deduped.realigned.recalibrated.mutect2.filtered.vcf.gz")
		String plasmabasename = basename("~{plasmabam}", ".filter.deduped.realigned.recalibrated.bam")
		String modules = "mrdetect/1.0 bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0 tabix"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
		String tumorVCFfilter = "FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'" 
		String tumorVAF = "0.1"
		String pickle = "$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/trained_SVM.pkl"
		String blacklist = "$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/blacklist.txt.gz"
		String genome = "$HG38_ROOT/hg38_random.fa"
		String difficultRegions = "--regions-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"
		String filterAndDetectScript = "$MRDETECT_ROOT/bin/filterAndDetect"
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		tumorvcf: "name of tumor sample in vcf"
		tumorbasename: "Base name for tumor"
		plasmabasename: "Base name for plasma"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		tumorVCFfilter: "set of filter calls to incl. in tumor VCF (any line with these flags will be included"
		tumorVAF: "Variant Allele Frequency for tumor VCF"
		pickle: "trained pickle for detecting real tumor reads"
		blacklist: "list of sites to exclude from analysis, gzipped"
		genome: "Path to loaded genome .fa"
		difficultRegions: "Path to .bed excluding difficult regions, string must include the flag --regions-file "
		filterAndDetectScript: "location of filter and detect script"
	}

	command <<<
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

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File? snvDetectionFinalResult = "~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv"
		File snvDetectionReadsScored = "~{plasmabasename}_PLASMA_VS_TUMOR.svm.tsv"
	}

	meta {
		output_meta: {
			snvDetectionFinalResult: "Final result and call from SNV detection",
			snvDetectionReadsScored: "Reads with potential for tumor, with their scores"
		}
	}
}


task parseControls {
	input {
		File controlFileList
		String controlFileListLoc = "~{controlFileList}"
		Int jobMemory = 4
		Int timeout = 12
	}

	parameter_meta {
		controlFileList: "file with list of control files"
		controlFileListLoc: "location of file with list of control files, for python intake"
		jobMemory: "Memory for this task in GB"
		timeout: "Timeout in hours, needed to override imposed limits"
	}

	command <<<
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

	runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
	}

	output {
		Array[Array[String]] controlFiles = read_tsv(stdout()) 
	}
}


task snvDetectionSummary {
	input {
		File? sampleCalls
		Array[File] controlCalls
		String DetectionRScript = "$MRDETECT_SCRIPTS_ROOT/bin/pwg_test.R"
		String samplebasename = basename("~{sampleCalls}", ".PLASMA_VS_TUMOR_RESULT.csv")
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
		String modules = "mrdetect-scripts/1.0"
	}

	parameter_meta {
		sampleCalls: "file of detection rate call for sample"
		controlCalls: "array of file of detection rate calls for HBCs"
		samplebasename: "base name for files"
		DetectionRScript: "location of pwg_test.R"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		cat ~{sep=' ' controlCalls} | awk '$1 !~ "BAM" {print}' >~{samplebasename}.HBCs.txt

		Rscript --vanilla ~{DetectionRScript} -s ~{sampleCalls} -S ~{samplebasename} -c ~{samplebasename}.HBCs.txt

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File pWGS_svg = "~{samplebasename}.pWGS.svg"
		File HBC_calls = "~{samplebasename}.HBCs.txt"
	}

	meta {
		output_meta: {
			pWGS_svg : "JSON file of mrdetect results"
		}
	}
}


