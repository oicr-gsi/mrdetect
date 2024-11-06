version 1.0

struct genomeResources {
    String ref_fasta
    String filterVCF_modules
    String filterVCF_difficultRegions
	String detectSNVs_modules
}

struct controlResources {
    String parseControls_modules
    String parseControls_controlFileList
}

workflow mrdetect {
	input {
		File? plasmabam
		File? plasmabai
		String? plasmaSampleName 
		String tumorSampleName
		File tumorvcf
		File tumorvcfindex
		String reference
		String instrument
		Boolean full_analysis_mode = true
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		plasmaSampleName: "name for plasma sample (from bam)"
		tumorvcf: "tumor vcf file, bgzip"
		tumorvcfindex: "tumor vcf index file"
		tumorSampleName: "ID for WGS tumor sample, must match .vcf header"
		reference: "genome reference build. Only hg38 supported"
		instrument: "sequencing instrument used (Illumina NovaSeq X Plus or Illumina NovaSeq 6000)"
        full_analysis_mode: "Enable full analysis mode with this flag"
	}

	Map[String,genomeResources] resources = {
	"hg38": {
		"ref_fasta": "$HG38_ROOT/hg38_random.fa",
		"filterVCF_modules": "bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0",
		"filterVCF_difficultRegions": "$HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed",
		"detectSNVs_modules" : "mrdetect/2.0.0 pwgs-blocklist/hg38.1"
		}
	}

	Map[String,controlResources] controls = {
	"Illumina NovaSeq X Plus": {
		"parseControls_modules": "pwgs-hbc/2.0",
		"parseControls_controlFileList" : "PWGS_HBC_LIST"
		},
	"Illumina NovaSeq 6000": {
		"parseControls_modules": "pwgs-hbc/1.0",
		"parseControls_controlFileList" : "PWGS_HBC_LIST"
		}

	}

	call filterVCF {
		input:
		modules = resources[reference].filterVCF_modules,
		genome = resources[reference].ref_fasta,
		difficultRegions = resources[reference].filterVCF_difficultRegions,
		tumorvcf = tumorvcf,
		tumorvcfindex = tumorvcfindex,
		tumorSampleName = tumorSampleName
	}

	if(full_analysis_mode) {
		
		call parseControls {
			input:
			modules = controls[instrument].parseControls_modules,
			controlFileList = controls[instrument].parseControls_controlFileList
		}


		scatter (control in parseControls.controlFiles) {
			call pullReads as controlPullReads {
				input:
				plasmabam = control[0],
				plasmabai = control[1],
				tumorvcf = filterVCF.filteredvcf,
				modules = resources[reference].detectSNVs_modules
			}
			call calculateQualityScore as controlQualityScore {
				input:
				snvDetectionReads = controlPullReads.snvDetectionReads,
				modules = resources[reference].detectSNVs_modules
			}
			call detectSNVs as controlDetectSNVs {
				input:
				plasmaSampleName = basename(control[0], ".bam"),
				tumorvcf = filterVCF.filteredvcf,
				snvDetectionReadsScored = controlQualityScore.snvDetectionReadsScored,
				modules = resources[reference].detectSNVs_modules
			}
		}

		call pullReads as samplePullReads {
			input:
			plasmabam = plasmabam,
			plasmabai = plasmabai,
			tumorvcf = filterVCF.filteredvcf,
			modules = resources[reference].detectSNVs_modules
		}
		call calculateQualityScore as sampleQualityScore {
			input:
			snvDetectionReads = samplePullReads.snvDetectionReads,
			modules = resources[reference].detectSNVs_modules
		}
		call detectSNVs as sampleDetectSNVs {
			input:
			plasmaSampleName = plasmaSampleName,
			tumorvcf = filterVCF.filteredvcf,
			snvDetectionReadsScored = sampleQualityScore.snvDetectionReadsScored,
			modules = resources[reference].detectSNVs_modules
		}

		call snvDetectionSummary {
			input:
			controlCalls = select_all(controlDetectSNVs.snvDetectionFinalResult),
			sampleCalls = sampleDetectSNVs.snvDetectionFinalResult,
			snpcount = filterVCF.snpcount,
			vafFile = sampleDetectSNVs.snvDetectionVAF,
			plasmaSampleName = plasmaSampleName
		}
	}

	meta {
		author: "Felix Beaudry"
		email: "fbeaudry@oicr.on.ca"
		description: "Workflow for MRdetect, detection of Minimal Residual Disease from paired tumor-plasma sample"
		    dependencies: [
			{
			name: "mrdetect/1.0",
			url: "https://ctl.cornell.edu/technology/mrdetect-license-request/"
			},
			{
			name: "bcftools/1.9",
			url: "https://github.com/samtools/bcftools"
			}]
		    output_meta: {
		    pWGS_svg: {
			description: "pWGS svg",
			vidarr_label: "pWGS_svg"
		    },
		    snvDetectionResult: {
			description: "Result from SNV detection incl sample HBCs",
			vidarr_label: "snvDetectionResult"
		    },
		    final_call: {
			description: "Final file of mrdetect results",
			vidarr_label: "final_call"
		    },
		    snvDetectionVAF: {
			description: "VAF from SNV detection for sample",
			vidarr_label: "snvDetectionVAF"
		    },
		    snpcount: {
			description: "number of SNPs in vcf after filtering",
			vidarr_label: "snpcount"
		    }
		}
		}
	output {
		File? snvDetectionResult = snvDetectionSummary.all_calls
		File? pWGS_svg = snvDetectionSummary.pWGS_svg
		File snpcount = filterVCF.snpcount
		File? snvDetectionVAF = sampleDetectSNVs.snvDetectionVAF
		File? final_call = snvDetectionSummary.final_call
		File? filteredvcf = filterVCF.filteredvcf
	}
}

task filterVCF {
	input {
		File tumorvcf
		File tumorvcfindex
		String tumorSampleName
		String tumorVCFfilter = "FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"
		String tumorVAF = "0.1"
		String genome 
		String difficultRegions
		String modules
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		tumorvcf: "tumor vcf file, bgzip"
		tumorvcfindex: "tumor vcf index file"
		tumorSampleName: "ID for WGS tumor sample"
		tumorVCFfilter: "set of filter calls to exclude in tumor VCF (any line with these flags will be excluded"
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

		$BCFTOOLS_ROOT/bin/bcftools view -s ~{tumorSampleName} --regions-file ~{difficultRegions} ~{tumorvcf} |\
		$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome} |\
		$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='snps'" |\
		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{tumorVCFfilter}" |\
		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{tumorVAF}" > ~{tumorSampleName}.SNP.vcf

		awk '$1 !~ "#" {print}' ~{tumorSampleName}.SNP.vcf | wc -l > ~{tumorSampleName}.SNP.count.txt

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File filteredvcf = "~{tumorSampleName}.SNP.vcf"
		File snpcount = "~{tumorSampleName}.SNP.count.txt"
	}

	meta {
		output_meta: {
			filteredvcf: "Filtered vcf",
			snpcount: "Number of SNPs after filtering"
		}
	}
}

task pullReads {
	input {
		File? plasmabam
		File? plasmabai
		File tumorvcf
		String modules
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 20
		String pullreadsScript = "$MRDETECT_ROOT/bin/pull_reads"
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		tumorvcf: "filtered tumor vcf file"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		pullreadsScript: "pull_reads.py executable"
	}

	command <<<
		set -euo pipefail

		~{pullreadsScript} \
			--bam ~{plasmabam} \
			--vcf ~{tumorvcf} \
			--out PLASMA_VS_TUMOR.tsv

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File? snvDetectionReads = "PLASMA_VS_TUMOR.tsv" 
	}

	meta {
		output_meta: {
			snvDetectionReads: "Reads with potential for tumor"
		}
	}
}

task calculateQualityScore {
	input {
		File? snvDetectionReads
		String modules
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 20
		String pickle = "$MRDETECT_ROOT/bin/MRDetectSNV/trained_SVM.pkl"
		String qualityscoreScript = "$MRDETECT_ROOT/bin/quality_score"
	}

	parameter_meta {
		snvDetectionReads: "Reads with potential for tumor"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		pickle: "trained pickle for detecting real tumor reads"
		qualityscoreScript: "quality_score.py executable"
	}

	command <<<
		set -euo pipefail

		~{qualityscoreScript} \
			--pickle-name ~{pickle} \
			--detections ~{snvDetectionReads} \
			--output_file PLASMA_VS_TUMOR.svm.tsv

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File? snvDetectionReadsScored = "PLASMA_VS_TUMOR.svm.tsv" 
	}

	meta {
		output_meta: {
			snvDetectionReadsScored: "Reads with potential for tumor, with their scores",
		}
	}
}

task detectSNVs {
	input {
		File tumorvcf
		String? plasmaSampleName 
		String tumorSampleName = basename(tumorvcf, ".vcf")
		File? snvDetectionReadsScored
		String modules
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 20
		String blocklist = "$PWGS_BLOCKLIST_ROOT/blocklist.vcf.gz"
		String filterAndDetectScript = "$MRDETECT_ROOT/bin/filterAndDetect"
	}

	parameter_meta {
		tumorvcf: "filtered tumor vcf file"
		plasmaSampleName: "name for plasma sample (from bam)"
		tumorSampleName: "name for tumour sample (from vcf)"
		blocklist: "list of sites to exclude from analysis, gzipped"
		snvDetectionReadsScored: "Reads with potential for tumor, with their scores"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		filterAndDetectScript: "filterAndDetect.py executable"
	}

	command <<<
		set -euo pipefail

		~{filterAndDetectScript} \
			--vcfid ~{tumorSampleName} --bamid ~{plasmaSampleName} \
			--svm ~{snvDetectionReadsScored} \
			--vcf ~{tumorvcf} \
			--output ./ \
			--blocklist ~{blocklist} \
	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File? snvDetectionFinalResult = "~{plasmaSampleName}.mrdetect.results.csv" #CHECK?
		File? snvDetectionVAF = "~{plasmaSampleName}.mrdetect.vaf.txt"
	}

	meta {
		output_meta: {
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
		String modules
	}

	parameter_meta {
		controlFileList: "file with list of control files"
		jobMemory: "Memory for this task in GB"
		timeout: "Timeout in hours, needed to override imposed limits"
	}

	command <<<
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

	runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
		modules: "~{modules}"
        }

	output {
		Array[Array[File]] controlFiles = read_tsv(stdout())
	}
}

task snvDetectionSummary {
	input {
		File? sampleCalls
		Array[File] controlCalls
		File snpcount
		File? vafFile
		String pvalue = "0.00001"
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
		String modules = "mrdetect/1.1.1"
		String pwgtestscript = "$MRDETECT_ROOT/bin/pwg_test"
		String? plasmaSampleName
	}

	parameter_meta {
		sampleCalls: "file of detection rate call for sample"
		controlCalls: "array of file of detection rate calls for HBCs"
		snpcount: "count of candidate SNPs"
		vafFile: "vaf from primary plasma"
		pvalue: "p-value for HBC error rate"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		pwgtestscript: "executable of pwg_test.R"
		plasmaSampleName: "name for plasma sample (from bam)"
	}

	command <<<
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

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File pWGS_svg = "~{plasmaSampleName}.pWGS.svg"
		File all_calls = "~{plasmaSampleName}.HBCs.csv"
		File final_call = "~{plasmaSampleName}.mrdetect.txt"
	}

	meta {
		output_meta: {
			pWGS_svg : "SVG plot of mrdetect results",
			all_calls : "HBC and sample mrdetect results",
			final_call : "final result"
		}
	}
}
