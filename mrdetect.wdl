version 1.0

workflow mrdetect {
	input {
		File plasmabam 
		File plasmabai 
		File normalbam 
		File normalbai 
		File tumorvcf
		String plasmabasename = basename("~{plasmabam}", ".bam")
		File segFile
		String plasmaSampleName
		String normalSampleName
		Int window = 500 
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		normalbam: "normal-tissue input .bam file"
		normalbai: "normal-tissue input .bai file"
		tumorvcf: "tumor vcf file"
		plasmabasename: "Base name for plasma"
		segFile: "segments file, eg from sequenza "
		plasmaSampleName: "plasma sample name"
		normalSampleName: "normal/reference sample name"
		window: "window size for scanning within intervals"
	}

	call detectSNVs {
		input: 
		plasmabam = plasmabam, 
		plasmabai = plasmabai, 
		tumorvcf = tumorvcf
	}

	call segTObed {
		input:
		segFile = segFile
	}

	call expandRegions { 
		input: 
		bedIntervals = segTObed.bedFile
	}

	scatter ( r in expandRegions.regions ) {
		call DepthOfCoverage as plasmaDepth {
			input:
			inputbam = plasmabam,
			inputbai = plasmabai,
			region = r[0],
			interval = r[1],
			CNVcall = r[2],
			sampleName = plasmaSampleName
  		}
  		call DepthOfCoverage as normalDepth {
			input:
			inputbam = normalbam,
			inputbai = normalbai,
			region = r[0],
			interval = r[1],
			CNVcall = r[2],
			sampleName = normalSampleName
  		}
	}

	call calculateMedians as plasmaMedians {
		input:
		depthOfCoverages = select_all(plasmaDepth.DepthOfCoverageOut),
		bedIntervals = segTObed.bedFile,
		sampleName = plasmaSampleName,
		window = window
	}

	call calculateMedians as normalMedians {
		input:
		depthOfCoverages = select_all(normalDepth.DepthOfCoverageOut),
		bedIntervals = segTObed.bedFile,
		sampleName = normalSampleName,
		window = window
	}

	call detectCNAs {
		input:
		normalCNVMedians = normalMedians.delsdupsMedian,
		normalNEUMedians = normalMedians.neuMedian,
		plasmaCNVMedians = plasmaMedians.delsdupsMedian,
		plasmaNEUMedians = plasmaMedians.neuMedian,
		normalSampleName = normalSampleName,
		plasmaSampleName = plasmaSampleName,
		window = window
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
			delsdupsZscore: "Final result and call from CNA detection"
		}
	}
	output {
		File snvDetectionFinalResult = "~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv"
		File delsdupsZscore = "~{plasmaSampleName}.~{normalSampleName}.win~{window}.robustZscore.delsdups.summary.txt"
	}
}

task detectSNVs {
	input {
		File plasmabam 
		File plasmabai 
		File tumorvcf
		String tumorbasename = basename("~{tumorvcf}", ".vcf.gz")
		String plasmabasename = basename("~{plasmabam}", ".bam")
		String modules = "mrdetect/1.0 bcftools/1.9"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
		String normalVCFfilter = "PASS,clustered_events"
		String normalVAF = "0.05"
		String pickle = "$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/trained_SVM.pkl"
		String blacklist = "$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/blacklist.txt.gz"
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
		normalVCFfilter: "set of filter calls to incl. in normal VCF (any line with these flags will be included"
		normalVAF: "Variant Allele Frequency for normal VCF"
		pickle: "trained pickle for detecting real tumor reads"
		blacklist: "list of sites to exclude from analysis, gzipped"
	}

	command <<<
		set -euo pipefail

		$BCFTOOLS_ROOT/bin/bcftools view \
			-f '~{normalVCFfilter}' \
			-v snps ~{tumorvcf}  | \
		$BCFTOOLS_ROOT/bin/bcftools filter \
			-i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{normalVAF}" >~{tumorbasename}.SNP.actuallyFiltered.vcf

		$MRDETECT_ROOT/bin/pull_reads \
			--bam ~{plasmabam} \
			--vcf ~{tumorbasename}.SNP.actuallyFiltered.vcf \
			--out ~{plasmabasename}_PLASMA_VS_TUMOR

		$MRDETECT_ROOT/bin/quality_score \
			--pickle-name ~{pickle} \
			--detections ~{plasmabasename}_PLASMA_VS_TUMOR.tsv \
			--output_file ~{plasmabasename}_PLASMA_VS_TUMOR.svm.tsv

		cp ~{blacklist} ./blacklist.txt.gz

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
			snvDetectionReadsScored: "Reads with potential for tumor, with their scores"
		}
	}
}

task segTObed {
	input {
		File segFile
		String segTObedrScript
		String segFileLoc = "~{segFile}"
		String basename = basename("~{segFile}", ".seg")
		Int jobMemory = 4
		Int timeout = 12
		String modules = "rstats/4.0"
	}
	parameter_meta{
		segFile: "segments file, eg from sequenza"
		segFileLoc: "segments file location as string, for R intake"
		basename: "Base name for segment file"
		jobMemory: "Memory allocated for this job (GB)"
		timeout: "Hours before task timeout"
		modules: "Required environment modules"
		segTObedrScript: "seg to bed R script location"
	}
	command <<<
		set -euo pipefail

		Rscript --vanilla ~{segTObedrScript} ~{segFileLoc} ~{basename}
	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
	}

	output {
		File bedFile = "~{basename}.seg.bed"
	}
}

task expandRegions {
	input {
		File bedIntervals
		String bedIntervalsLoc = "~{bedIntervals}"
		Int jobMemory = 4
		Int timeout = 12
	}

	parameter_meta {
		bedIntervals: "bed file with intervals"
		bedIntervalsLoc: "bed file location as string, for python intake"
		jobMemory: "Memory for this task in GB"
		timeout: "Timeout in hours, needed to override imposed limits"
	}

	command <<<
		python <<CODE
		import os, re

		with open("~{bedIntervalsLoc}") as f:
			for line in f:
				if re.match("#",line):
					pass
				else:
					line = line.rstrip()
					tmp = line.split("\t")
					r = tmp[0] + ":" + tmp[1] + "-" + tmp[2] + "\t" + tmp[0] + "_" + tmp[1] + "_" + tmp[2] + "\t" + tmp[3]
					print(r)
		f.close()
		CODE
	>>>

	runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
	}

	output {
		Array[Array[String]] regions = read_tsv(stdout()) 
	}
}

task DepthOfCoverage {
	input {
		File inputbam
		File inputbai
		String sampleName
		String region
		String interval
		String CNVcall
		String ref = "$HG38_ROOT/hg38_random.fa"
		String modules = "mrdetect/1.0 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		inputbam: " input .bam file"
		inputbai: " input .bam file"
		sampleName: "sample Name"
		region: "Region in form chrX:12000-12500"
		interval: "Region in form chrX_12000_12500"
		CNVcall: "DUP, DEL or NEU"
		modules: "Required environment modules"
		ref: "Genome Reference"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		java -Xmx24576m -jar $GATK_ROOT/GenomeAnalysisTK.jar \
			-T DepthOfCoverage \
			--includeDeletions \
			--read_filter BadCigar \
			--summaryCoverageThreshold 100 \
			--intervals ~{region} \
			--interval_merging OVERLAPPING_ONLY \
			-R ~{ref} \
			-I ~{inputbam} \
			-o ~{sampleName}.DepthOfCoverage

		mv ~{sampleName}.DepthOfCoverage ~{interval}.~{CNVcall}.~{sampleName}.DepthOfCoverage 

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File? DepthOfCoverageOut = "~{interval}.~{CNVcall}.~{sampleName}.DepthOfCoverage"
	}

	meta {
		output_meta: {
			DepthOfCoverageOut: "Depth of Coverage output for interval"
		}
	}
}

task calculateMedians {
	input {
		Array[File] depthOfCoverages
		String sampleName
		Int window
		File bedIntervals
		String modules = "mrdetect/1.0 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		depthOfCoverages: "depth of coverage in plasma for each segment"
		sampleName: "sample Name"
		window: "window size for scanning within intervals"
		bedIntervals: "bed file with intervals"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		mkdir segments
		ln -s ~{sep=' ' depthOfCoverages} segments/

		$MRDETECT_ROOT/bin/printmedian \
		-i segments/ -o segments/ \
		-s ~{sampleName} \
		-f ~{bedIntervals} \
		-w ~{window} -t CNV -z

		$MRDETECT_ROOT/bin/printmedian \
		-i segments/ -o segments/ \
		-s ~{sampleName} \
		-f ~{bedIntervals} \
		-w ~{window} -t NEU -z

		mv segments/~{sampleName}.win~{window}.delsdups.perwindow.median ~{sampleName}.win~{window}.delsdups.perwindow.median 
		mv segments/~{sampleName}.win~{window}.neu.perwindow.median ~{sampleName}.win~{window}.neu.perwindow.median 

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File delsdupsMedian = "~{sampleName}.win~{window}.delsdups.perwindow.median"
		File neuMedian = "~{sampleName}.win~{window}.neu.perwindow.median"
	}

	meta {
		output_meta: {
			delsdupsMedian: "Median depth of coverage for deletions and duplicates",
			neuMedian: "Median depth of coverage for neutral segments"
		}
	}
}

task detectCNAs {
	input {
		File plasmaCNVMedians
		File normalCNVMedians
		File plasmaNEUMedians
		File normalNEUMedians
		String normalSampleName
		String plasmaSampleName
		Int tumor_coverage 
		Int normal_coverage
		Int window
		Float? median_thresh = 1.5
		String? interval_name = "interval"
		String modules = "mrdetect/1.0 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		plasmaCNVMedians: "Median depth of coverage for deletions and duplicates in plasma"
		normalCNVMedians: "Median depth of coverage for deletions and duplicates in normal"
		plasmaNEUMedians: "Median depth of coverage for neutral segments in plasma"
		normalNEUMedians: "Median depth of coverage for neutral segments in normal"
		normalSampleName: "normal sample Name"
		plasmaSampleName: "plasma sample Name"
		tumor_coverage: "mean coverage in tumor, for normalization"
		normal_coverage: "mean coverage in normal, for normalization"
		window: "window size for scanning within intervals"
		median_thresh: "median threshold, removed bins with extreme coverage values (>abs(1.5*median))"
		interval_name: "name for interval"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		mkdir input input/~{normalSampleName} input/~{normalSampleName}/~{interval_name} input/~{normalSampleName}/~{interval_name}/segments input/~{plasmaSampleName} input/~{plasmaSampleName}/~{interval_name} input/~{plasmaSampleName}/~{interval_name}/segments output

		ln -s ~{normalCNVMedians} ~{normalNEUMedians} input/~{normalSampleName}/~{interval_name}/segments
		ln -s ~{plasmaCNVMedians} ~{plasmaNEUMedians} input/~{plasmaSampleName}/~{interval_name}/segments

		$MRDETECT_ROOT/bin/parsemedian_robustZscore -i input -t ~{plasmaSampleName} -n ~{normalSampleName} -in ~{interval_name} -o output -w ~{window} -mt ~{median_thresh} -tc ~{tumor_coverage}  -nc ~{normal_coverage}

		cp output/~{plasmaSampleName}.~{normalSampleName}.win~{window}.robustZscore.delsdups.summary.txt ~{plasmaSampleName}.~{normalSampleName}.win~{window}.robustZscore.delsdups.summary.txt
	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File delsdupsZscore = "~{plasmaSampleName}.~{normalSampleName}.win~{window}.robustZscore.delsdups.summary.txt"
	}

	meta {
		output_meta: {
			delsdupsZscore: "Z-score for detection of deletions and duplicates"
		}
	}
}
