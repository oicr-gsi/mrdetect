version 1.0

workflow mrdetect {
	input {
		File plasmabam 
		File plasmabai
		String plasmabasename = basename("~{plasmabam}", ".filter.deduped.realigned.recalibrated.bam")
		File tumorvcf
		File controlFileList
		Boolean CNVtest
		File? segFile
		File? referencebam 
		File? referencebai
		String? referencebasename = basename("~{referencebam}", ".filter.deduped.realigned.recalibrated.bam")
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		referencebam: "reference (or normal) input .bam file"
		referencebai: "reference (or normal) input .bai file"
		tumorvcf: "tumor vcf file"
		plasmabasename: "Base name for plasma"
		CNVtest: "should CNV analysis be run?"
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

	if ( CNVtest == true && defined(referencebam) && defined(referencebai) && defined(segFile) && defined(referencebasename) ) {

		call calcCoverage as plasmaCoverage {
			input:
			inputbam = plasmabam
		}

		call calcCoverage as referenceCoverage {
			input:
			inputbam = referencebam
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
				CNVcall = r[2]
	  		}
	  		call DepthOfCoverage as referenceDepth {
				input:
				inputbam = referencebam,
				inputbai = referencebai,
				region = r[0],
				interval = r[1],
				CNVcall = r[2]
	  		}
		}

		call calculateMedians as plasmaMedians {
			input:
			depthOfCoverages = select_all(plasmaDepth.DepthOfCoverageOut),
			bedIntervals = segTObed.bedFile,
			basename = plasmabasename
		}

		call calculateMedians as referenceMedians {
			input:
			depthOfCoverages = select_all(referenceDepth.DepthOfCoverageOut),
			bedIntervals = segTObed.bedFile,
			basename = referencebasename
		}

		call detectCNAs {
			input:
			referenceCNVMedians = referenceMedians.delsdupsMedian,
			referenceNEUMedians = referenceMedians.neuMedian,
			plasmaCNVMedians = plasmaMedians.delsdupsMedian,
			plasmaNEUMedians = plasmaMedians.neuMedian,
			plasma_coverage = plasmaCoverage.meanCoverage,
			reference_coverage = referenceCoverage.meanCoverage,
			referencebasename = referencebasename,
			plasmabasename = plasmabasename
		}
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
			delsdupsZscore: "Final result and call from CNA detection",
			pWGS_svg: "pWGS svg"
		}
	}
	output {
		File snvDetectionFinalResult = "~{plasmabasename}_PLASMA_VS_TUMOR_RESULT.csv"
		File? delsdupsZscore = "~{plasmabasename}.robustZscore.delsdups.summary.txt"
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
		String modules = "mrdetect/1.0 bcftools/1.9 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
		String tumorVCFfilter = "FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'" 
		String tumorVAF = "0.1"
		String pickle = "$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/trained_SVM.pkl"
		String blacklist = "$MRDETECT_ROOT/MRDetect-master/MRDetectSNV/blacklist.txt.gz"
		String genome = "$HG38_ROOT/hg38_random.fa"

	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bam file"
		tumorvcf: "tumor vcf file"
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
	}

	command <<<
		set -euo pipefail

		$BCFTOOLS_ROOT/bin/bcftools view -s ~{tumorbasename} ~{tumorvcf} |\
		$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome} |\
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

		$MRDETECT_ROOT/bin/filterAndDetect \
			~{tumorbasename}.SNP.vcf \
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
		String? samplebasename = basename("~{sampleCalls}", ".PLASMA_VS_TUMOR_RESULT.csv")
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
		String modules = "mrdetect-scripts/1.0"
	}

	parameter_meta {
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

task calcCoverage {
	input {
		File? inputbam
		String modules = "samtools"
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		samtools depth -a ~{inputbam} |  awk '{sum+=$3} END {print sum/NR}'

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		Int meanCoverage = stdout()
	}

	meta {
		output_meta: {
			JSONout : "JSON file of sigtools and CHORD signatures"
		}
	}
}

task segTObed {
	input {
		File? segFile
		String segTObedrScript = "/.mounts/labs/CGI/scratch/fbeaudry/mrdetect/segTObed.R"
		String segFileLoc = "~{segFile}"
		String basename = basename("~{segFile}", ".seg")
		Int jobMemory = 4
		Int timeout = 12
		String modules = "rstats/4.0"
	}
	parameter_meta{
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
		File? inputbam
		File? inputbai
		String basename = basename("~{inputbam}", ".filter.deduped.realigned.recalibrated.bam")
		String region
		String interval
		String CNVcall
		String genome = "$HG38_ROOT/hg38_random.fa"
		String modules = "mrdetect/1.0 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		inputbam: " input .bam file"
		inputbai: " input .bam file"
		basename: "sample Name"
		region: "Region in form chrX:12000-12500"
		interval: "Region in form chrX_12000_12500"
		CNVcall: "DUP, DEL or NEU"
		modules: "Required environment modules"
		genome: "Genome Reference"
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
			-R ~{genome} \
			-I ~{inputbam} \
			-o ~{basename}.DepthOfCoverage

		mv ~{basename}.DepthOfCoverage ~{interval}.~{CNVcall}.~{basename}.DepthOfCoverage 

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File? DepthOfCoverageOut = "~{interval}.~{CNVcall}.~{basename}.DepthOfCoverage"
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
		String? basename
		Int window = 500
		File bedIntervals
		String modules = "mrdetect/1.0 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		depthOfCoverages: "depth of coverage in plasma for each segment"
		basename: "sample Name"
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
		-s ~{basename} \
		-f ~{bedIntervals} \
		-w ~{window} -t CNV -z

		$MRDETECT_ROOT/bin/printmedian \
		-i segments/ -o segments/ \
		-s ~{basename} \
		-f ~{bedIntervals} \
		-w ~{window} -t NEU -z

		mv segments/~{basename}.delsdups.perwindow.median ~{basename}.delsdups.perwindow.median 
		mv segments/~{basename}.neu.perwindow.median ~{basename}.neu.perwindow.median 

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File delsdupsMedian = "~{basename}.delsdups.perwindow.median"
		File neuMedian = "~{basename}.neu.perwindow.median"
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
		File referenceCNVMedians
		File plasmaNEUMedians
		File referenceNEUMedians
		String? referencebasename
		String plasmabasename
		Int plasma_coverage 
		Int reference_coverage
		Int window  = 500
		Float? median_thresh = 1.5
		String? interval_name = "interval"
		String modules = "mrdetect/1.0 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		plasmaCNVMedians: "Median depth of coverage for deletions and duplicates in plasma"
		referenceCNVMedians: "Median depth of coverage for deletions and duplicates in reference (normal)"
		plasmaNEUMedians: "Median depth of coverage for neutral segments in plasma"
		referenceNEUMedians: "Median depth of coverage for neutral segments in reference (normal)"
		referencebasename: "reference sample Name"
		plasmabasename: "plasma sample Name"
		plasma_coverage: "mean coverage in plasma, for normalization"
		reference_coverage: "mean coverage in reference, for normalization"
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

		mkdir input input/~{referencebasename} input/~{referencebasename}/~{interval_name} input/~{referencebasename}/~{interval_name}/segments input/~{plasmabasename} input/~{plasmabasename}/~{interval_name} input/~{plasmabasename}/~{interval_name}/segments output

		ln -s ~{referenceCNVMedians} ~{referenceNEUMedians} input/~{plasmabasename}/~{interval_name}/segments
		ln -s ~{plasmaCNVMedians} ~{plasmaNEUMedians} input/~{plasmabasename}/~{interval_name}/segments

		$MRDETECT_ROOT/bin/parsemedian_robustZscore -i input -t ~{plasmabasename} -n ~{referencebasename} -in ~{interval_name} -o output -w ~{window} -mt ~{median_thresh} -tc ~{plasma_coverage} -nc ~{reference_coverage}

		cp output/~{plasmabasename}.~{referencebasename}.win~{window}.robustZscore.delsdups.summary.txt ~{plasmabasename}.robustZscore.delsdups.summary.txt
	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File delsdupsZscore = "~{plasmabasename}.robustZscore.delsdups.summary.txt"
	}

	meta {
		output_meta: {
			delsdupsZscore: "Z-score for detection of deletions and duplicates"
		}
	}
}
