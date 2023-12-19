#!/usr/bin/python3

# create a new output file
outfile = open('../../refs/config.json', 'w')

# get all file names
allSamples = list()
read = ["R1", "R2"]
numSamples = 0

with open('../../refs/sample_file_list.txt', 'r') as infile:
    for line in infile:
        numSamples += 1
        sample = line.strip()
        allSamples.append(sample.replace("_R1.fastq.gz", ""))

# create header and write to outfile
header = '''{{
    "DIRECTORIES",
    "rawReads" : "rawReads/",
    "rawQC" : "rawQC/",
    "trimmedReads" : "trimmedReads/",
    "trimmedQC" : "trimmedQC/",
    "starAligned" : "starAligned/",
    "featureCounts" : "featureCounts/",
    "genomeDir" : "refs/starGenomeDir/",

    "FILES",
    "Mmusculus.gtf" : "/research/labs/neurology/fryer/projects/references/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf",
    "Mmusculus.fa" : "/research/labs/neurology/fryer/projects/references/mouse/refdata-gex-mm10-2020-A/fasta/genome.fa",

    "SAMPLE INFORMATION",
    "allSamples": {0},

    "CLUSTER INFORMATION",
    "threads" : "10",
'''
outfile.write(header.format(allSamples))


# config formatting
counter = 0
with open('../../refs/sample_file_list.txt', 'r') as infile:
    for line in infile:
        counter += 1

        # store filename
        sample = line.strip()
        read1 = sample.replace(".fastq.gz", "")
        read2 = sample.replace("_R1.fastq.gz", "_R2")
        baseName = sample.replace("_R1.fastq.gz", "")

        # break down fastq file info
        # @A00127:312:HVNLJDSXY:2:1101:2211:1000
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>

        out = '''
    "{0}":{{
        "read1": "{1}",
        "read2": "{2}",
        '''
        outfile.write(out.format(baseName, read1, read2))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()

