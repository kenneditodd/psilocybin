configfile: "refs/config.json"


# DIRECTORY VARIABLES
#-----------------------------------------------------------------------------------------
rawReadsDir = config["rawReads"]
trimmedReadsDir = config["trimmedReads"]
rawQCDir = config["rawQC"]
trimmedQCDir = config["trimmedQC"]
starDir = config["starAligned"]
countsDir = config["featureCounts"]


# RULE ALL
#-----------------------------------------------------------------------------------------
rule all:
	input:
	  expand(config["genomeDir"]),
		expand(trimmedReadsDir + "{sample}_trim_R1.fastq.gz", sample = config["allSamples"]),
		expand(trimmedReadsDir + "{sample}_trim_R2.fastq.gz", sample = config["allSamples"]),
		expand(countsDir + "{sample}_gene.counts", sample = config["allSamples"]),
		expand(countsDir + "{sample}_exon.counts", sample = config["allSamples"])


# INDEX GENOME
#-----------------------------------------------------------------------------------------
rule index_genome:
	input:
		fa = config["Mmusculus.fa"],
		gtf = config["Mmusculus.gtf"]
	output:
		starIndex = directory(config["genomeDir"]),
	params:
		genomeDir = config["genomeDir"],
		threads = config["threads"]
	shell:
		"""
		STAR --runMode genomeGenerate --runThreadN {params.threads} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeDir {params.genomeDir}
		"""



# TRIM BBDUK
#-----------------------------------------------------------------------------------------
rule trim_bbduk:
	input:
		R1 = lambda wildcards: rawReadsDir + config[wildcards.sample]["read1"] + ".fastq.gz",
		R2 = lambda wildcards: rawReadsDir + config[wildcards.sample]["read2"] + ".fastq.gz"
	output:
		trimR1 = trimmedReadsDir + "{sample}_trim_R1.fastq.gz",
		trimR2 = trimmedReadsDir + "{sample}_trim_R2.fastq.gz"
	params:
	  threads = config["threads"]
	shell:
		"""
		bbduk.sh -Xmx3g in1={input.R1} in2={input.R2} out1={output.trimR1} out2={output.trimR2} ref=refs/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo threads={params.threads} trimpolyg=1 trimpolya=1
		"""


# ALIGN READS
#-----------------------------------------------------------------------------------------
rule align_reads:
	input:
		trimR1 = trimmedReadsDir + "{sample}_trim_R1.fastq.gz",
		trimR2 = trimmedReadsDir + "{sample}_trim_R2.fastq.gz",
		genomeDir = config["genomeDir"]
	output:
		aligned = (starDir + "{sample}.Aligned.sortedByCoord.out.bam")
	params:
		prefix = (starDir + "{sample}."),
		threads = config["threads"]
	shell:
		"""
		STAR --genomeDir {input.genomeDir} --runThreadN {params.threads} --readFilesCommand zcat --limitBAMsortRAM 31000000000 --readFilesIn {input.trimR1} {input.trimR2} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate
		"""


# FEATURE COUNTS
#-----------------------------------------------------------------------------------------	
rule gene_count:
	input:
		bam = starDir + "{sample}.Aligned.sortedByCoord.out.bam",
		gtf = config["Mmusculus.gtf"]
	output:
		feature = countsDir + "{sample}_gene.counts"
	params:
	  threads = config["threads"]
	shell:
		"""
		featureCounts -p --primary -t gene -T {params.threads} -s 2 -a {input.gtf} -o {output.feature} {input.bam}
		
		# KEY
		# -p specify that input data contains paired-end reads
		# --primary count primary alignments only, primary alignments are identified using bit 0x100 in SAM/BAM FLAG field
		# -t specify feature type in annotation, exon by default
		# -T number of the threads, 1 by default
		# -s specify strandedness, 0 = unstranded, 1 = stranded, 2 = reversely stranded
		# -a name of an annotation file. GTF/GFF format by default
		# -o name of output file including read counts
		"""

rule exon_count:
	input:
		bam = starDir + "{sample}.Aligned.sortedByCoord.out.bam",
		gtf = config["Mmusculus.gtf"]
	output:
		feature = countsDir + "{sample}_exon.counts"
	params:
	  threads = config["threads"]
	shell:
		"""
		featureCounts -p --primary -t exon -T {params.threads} -s 2 -a {input.gtf} -o {output.feature} {input.bam}
		
		# KEY
		# -p specify that input data contains paired-end reads
		# --primary count primary alignments only, primary alignments are identified using bit 0x100 in SAM/BAM FLAG field
		# -t specify feature type in annotation, exon by default
		# -T number of the threads, 1 by default
		# -s specify strandedness, 0 = unstranded, 1 = stranded, 2 = reversely stranded
		# -a name of an annotation file. GTF/GFF format by default
		# -o name of output file including read counts
		"""

