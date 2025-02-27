configfile: "refs/config.json"


# CONFIG VARIABLES
#--------------------------------------------------------------------------------
trimmedReadsDir = config["DIRECTORIES"]["trimmedReads"]
rawQCDir = config["DIRECTORIES"]["rawQC"]
trimmedQCDir = config["DIRECTORIES"]["trimmedQC"]
starDir = config["DIRECTORIES"]["starAligned"]
countsDir = config["DIRECTORIES"]["featureCounts"]


# ENVIRONMENT VARIABLES
#-------------------------------------------------------------------------------
import os
from dotenv import load_dotenv

# Load the .env file
load_dotenv("refs/.env")

# Get paths from environment variables
gtf_file = os.getenv("MMUSCULUS_GTF")
fa_file = os.getenv("MMUSCULUS_FA")
rawReadsDir = os.getenv("RAW_READS_DIR")


# RULE ALL
#-------------------------------------------------------------------------------
rule all:
    input:
        config["DIRECTORIES"]["genomeDir"] + "Genome",
        expand(trimmedReadsDir + "{sample}_trim_R1.fastq.gz", sample=config["SAMPLE_INFORMATION"]["allSamples"]),
        expand(trimmedReadsDir + "{sample}_trim_R2.fastq.gz", sample=config["SAMPLE_INFORMATION"]["allSamples"]),
        expand(rawQCDir + "{sample}_R1_001_fastqc.zip", sample=config["SAMPLE_INFORMATION"]["allSamples"]),
        expand(trimmedQCDir + "{sample}_trim_R1_fastqc.zip", sample=config["SAMPLE_INFORMATION"]["allSamples"]),
        rawQCDir + "multiqc_report.html",
        trimmedQCDir + "multiqc_report.html",
        expand(starDir + "{sample}.Aligned.sortedByCoord.out.bam", sample=config["SAMPLE_INFORMATION"]["allSamples"]),
        expand(countsDir + "{sample}_{feature_type}.counts", sample=config["SAMPLE_INFORMATION"]["allSamples"], feature_type=["gene", "exon"])


# INDEX GENOME
#-------------------------------------------------------------------------------
rule index_genome:
    input:
        fa = fa_file,
        gtf = gtf_file
    output:
        starIndex = config["DIRECTORIES"]["genomeDir"] + "Genome"
    params:
        genomeDir = config["DIRECTORIES"]["genomeDir"],
        threads = config["CLUSTER_INFORMATION"]["threads"]
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN {params.threads} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeDir {params.genomeDir}
        """


# RAW FASTQC
#-------------------------------------------------------------------------------
rule raw_fastqc:
    input:
        R1 = lambda wildcards: rawReadsDir + config["SAMPLES"][wildcards.sample]["read1"] + ".fastq.gz",
        R2 = lambda wildcards: rawReadsDir + config["SAMPLES"][wildcards.sample]["read2"] + ".fastq.gz"
    output:
        qc1 = rawQCDir + "{sample}_R1_001_fastqc.zip",
        qc2 = rawQCDir + "{sample}_R2_001_fastqc.zip"
    params:
        threads = config["CLUSTER_INFORMATION"]["threads"]
    shell:
        """
        fastqc {input.R1} {input.R2} -o {rawQCDir} --threads {params.threads}
        """


# TRIM BBDUK
#-------------------------------------------------------------------------------
rule trim_bbduk:
    input:
        R1 = lambda wildcards: rawReadsDir + config["SAMPLES"][wildcards.sample]["read1"] + ".fastq.gz",
        R2 = lambda wildcards: rawReadsDir + config["SAMPLES"][wildcards.sample]["read2"] + ".fastq.gz"
    output:
        trimR1 = trimmedReadsDir + "{sample}_trim_R1.fastq.gz",
        trimR2 = trimmedReadsDir + "{sample}_trim_R2.fastq.gz"
    params:
        threads = config["CLUSTER_INFORMATION"]["threads"]
    shell:
        """
        bbduk.sh -Xmx3g in1={input.R1} in2={input.R2} out1={output.trimR1} out2={output.trimR2} ref=refs/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo threads={params.threads} trimpolyg=1 trimpolya=1
        """


# TRIMMED FASTQC
#-------------------------------------------------------------------------------
rule trimmed_fastqc:
    input:
        R1 = trimmedReadsDir + "{sample}_trim_R1.fastq.gz",
        R2 = trimmedReadsDir + "{sample}_trim_R2.fastq.gz"
    output:
        qc1 = trimmedQCDir + "{sample}_trim_R1_fastqc.zip",
        qc2 = trimmedQCDir + "{sample}_trim_R2_fastqc.zip"
    params:
        threads = config["CLUSTER_INFORMATION"]["threads"]
    shell:
        """
        fastqc {input.R1} {input.R2} -o {trimmedQCDir} --threads {params.threads}
        """


# ALIGN READS
#-------------------------------------------------------------------------------
rule align_reads:
    input:
        trimR1 = trimmedReadsDir + "{sample}_trim_R1.fastq.gz",
        trimR2 = trimmedReadsDir + "{sample}_trim_R2.fastq.gz",
        genomeIndex = rules.index_genome.output.starIndex
    output:
        aligned = starDir + "{sample}.Aligned.sortedByCoord.out.bam"
    params:
        prefix = starDir + "{sample}.",
        threads = config["CLUSTER_INFORMATION"]["threads"],
        genomeDir = config["DIRECTORIES"]["genomeDir"]
    shell:
        """
        STAR --genomeDir {params.genomeDir} --runThreadN {params.threads} --readFilesCommand zcat --limitBAMsortRAM 31000000000 --readFilesIn {input.trimR1} {input.trimR2} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate
        """


# FEATURE COUNTS
#-------------------------------------------------------------------------------   
rule feature_count:
    input:
        bam = starDir + "{sample}.Aligned.sortedByCoord.out.bam",
        gtf = gtf_file
    output:
        counts = countsDir + "{sample}_{feature_type}.counts"
    params:
        feature_type = "{feature_type}",  # can be 'gene' or 'exon'
        threads = config["CLUSTER_INFORMATION"]["threads"]
    shell:
        """
        featureCounts -p --primary -t {params.feature_type} -T {params.threads} -s 2 -a {input.gtf} -o {output.counts} {input.bam}
        """


# MULTIQC FOR RAW FASTQC
#-------------------------------------------------------------------------------
rule raw_multiqc:
    input:
        expand(rawQCDir + "{sample}_R1_001_fastqc.zip", sample=config["SAMPLE_INFORMATION"]["allSamples"]),
        expand(rawQCDir + "{sample}_R2_001_fastqc.zip", sample=config["SAMPLE_INFORMATION"]["allSamples"])
    output:
        html = rawQCDir + "multiqc_report.html"
    params:
        outdir = rawQCDir
    shell:
        """
        multiqc --interactive {rawQCDir} -o {params.outdir}
        """


# MULTIQC FOR TRIMMED FASTQC
#-------------------------------------------------------------------------------
rule trimmed_multiqc:
    input:
        expand(trimmedQCDir + "{sample}_trim_R1_fastqc.zip", sample=config["SAMPLE_INFORMATION"]["allSamples"]),
        expand(trimmedQCDir + "{sample}_trim_R2_fastqc.zip", sample=config["SAMPLE_INFORMATION"]["allSamples"])
    output:
        html = trimmedQCDir + "multiqc_report.html"
    params:
        outdir = trimmedQCDir
    shell:
        """
        multiqc --interactive {trimmedQCDir} -o {params.outdir}
        """

