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
		plasmaSampleName: "plasma sample Name"
		normalSampleName: "normal/reference sample Name"
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
		bedIntervals = bedIntervals,
		sampleName = plasmaSampleName,
		window = window
	}

	call calculateMedians as normalMedians {
		input:
		depthOfCoverages = select_all(normalDepth.DepthOfCoverageOut),
		bedIntervals = bedIntervals,
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
			snvDetectionReadsScored: "Reads with potential for tumor, with their scores"
		}
	}
}

task segTObed {
	input {
		File segFile
		String segFileLoc = "~{segFile}"
		String basename = basename("~{segFile}", ".seg")
		Int jobMemory = 4
		Int timeout = 12
	}
	parameter_meta{

	}
	command <<<
		R <<CODE
		seg <- read.table('~{segFileLoc}')

		seg$type <- "NEU"
		seg$type[seg$seg.mean < -0.3] <- "DEL"
		seg$type[seg$seg.mean > 0.3] <- "DUP"

		seg <- seg[,c('chrom','loc.start','loc.end','type','seg.mean')]
		names(seg) <- c("#chr",  "start", "end", "type",  "log2")

		write.table(
		  seg,
		  file = "~{basename}.seg.bed",
		  append = F, quote = FALSE, sep = "\t", 
		  eol = "\n", na = "NA",dec = ".", row.names = FALSE, 
		  col.names = TRUE
		)
		CODE
	>>>

	runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
	}

	output {
		File bedFile = ~{basename}.seg.bed 
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
		String modules = "mrdetect/1.0 hg38/p12"
		String ref = "$HG38_ROOT/hg38_random.fa"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		inputbam: " input .bam file"
		inputbai: " input .bam file"
		modules: "Required environment modules"
		ref: "Genome Reference"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
		region: "Region in form chrX:12000-12500"
		interval: "Region in form chrX_12000_12500"
		CNVcall: "DUP, DEL or NEU"
		sampleName: "sample Name"
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
		File bedIntervals
		String modules = "mrdetect/1.0 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
		Int window
	}

	parameter_meta {
		window: "window size for scanning within intervals"
		sampleName: "sample Name"
		depthOfCoverages: "depth of coverage in plasma for each segment"
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
			delsdupsMedian: "Median depth of coverage for deletions and duplicates in plasma",
			neuMedian: "Median depth of coverage for neutral segments in plasma"
		}
	}
}

task detectCNAs {
	input {
		Float? median_thresh = 1.5
		Int tumor_coverage 
		Int normal_coverage
		Int window
		String? interval_name = "interval"
		File normalCNVMedians
		File normalNEUMedians
		File plasmaCNVMedians
		File plasmaNEUMedians
		String normalSampleName
		String plasmaSampleName
		String modules = "mrdetect/1.0 hg38/p12"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		window: "window size for scanning within intervals"
		normal_coverage: "mean coverage in normal, for normalization"
		tumor_coverage: "mean coverage in tumor, for normalization"
		median_thresh: "median threshold, removed bins with extreme coverage values (>abs(1.5*median))"
		interval_name: "name for interval"
		plasmaCNVMedians: "Median depth of coverage for deletions and duplicates in plasma"
		normalCNVMedians: "Median depth of coverage for deletions and duplicates in normal"
		plasmaNEUMedians: "Median depth of coverage for neutral segments in plasma"
		normalNEUMedians: "Median depth of coverage for neutral segments in normal"
		plasmaSampleName: "plasma sample Name"
		normalSampleName: "normal sample Name"
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
			delsdupsZscore: "Z-score for detection of deletions and duplicates",
		}
	}
}
