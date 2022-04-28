version 1.0

workflow mrdetect {
	input {
		File plasmabam 
		File plasmabai 
		File tumorvcf
		String plasmabasename = basename("~{plasmabam}", ".bam")
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		tumorvcf: "tumor vcf file"
		plasmabasename: "Base name for plasma"
	}

	call detectSNVs {
		input: 
		plasmabam = plasmabam, 
		plasmabai = plasmabai, 
		tumorvcf = tumorvcf
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
			},
			{
				name: "tabix/1.9",
				url: "http://www.htslib.org/doc/tabix.html"
			}
		]
		output_meta: {
			snvDetectionFinalResult: "Final result and call from SNV detection"
		}
	}
	output {
		File snvDetectionFinalResult = "~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv"
	}
}

task detectSNVs {
	input {
		File plasmabam 
		File plasmabai 
		File tumorvcf
		String tumorbasename = basename("~{tumorvcf}", ".vcf.gz")
		String plasmabasename = basename("~{plasmabam}", ".bam")
		String modules = "mrdetect/1.0 bcftools/1.9 tabix/1.9"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bam file"
		tumorvcf: "tumor vcf file"
		tumorbasename: "Base name for tumor"
		plasmabasename: "Base name for plasma"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		$BCFTOOLS_ROOT/bin/bcftools view \
			-f 'PASS,clustered_events' \
			-v snps ~{tumorvcf}  | \
		$BCFTOOLS_ROOT/bin/bcftools filter \
			-i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.05" >~{tumorbasename}.SNP.actuallyFiltered.vcf

		$MRDETECT_ROOT/bin/pull_reads \
			--bam ~{plasmabam} \
			--vcf ~{tumorbasename}.SNP.actuallyFiltered.vcf \
			--out ~{plasmabasename}_PLASMA_VS_TUMOR

		$MRDETECT_ROOT/bin/quality_score \
			--pickle-name $MRDETECT_ROOT/MRDetect-master/MRDetectSNV/trained_SVM.pkl \
			--detections ~{plasmabasename}_PLASMA_VS_TUMOR.tsv \
			--output_file ~{plasmabasename}_PLASMA_VS_TUMOR.svm.tsv

		cp $MRDETECT_ROOT/MRDetect-master/MRDetectSNV/blacklist.txt.gz ./blacklist.txt.gz

		$MRDETECT_ROOT/bin/filterAndDetect \
			~{tumorbasename}.SNP.actuallyFiltered.vcf \
			~{plasmabasename}_PLASMA_VS_TUMOR.svm.tsv \
			~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File snvDetectionFinalResult = "~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv"
		File snvDetectionReadsScored = "~{plasmabasename}_PLASMA_VS_TUMOR.svm.tsv"
	}

	meta {
		output_meta: {
			snvDetectionFinalResult: "Final result and call from SNV detection",
			snvDetectionReadsScored: "Reads with potential for tumor, with their scores",
		}
	}
}
