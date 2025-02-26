# Kennedi Todd
# September 17, 2024

# go to counts folder
cd ../featureCounts

# GENE COUNT MATRIX
#-------------------------------------------------------------------------------

# create list of gene featureCount count files
ls -1 | grep gene.counts$ > gene_count_files.txt

# store number of files in variable
n=$(wc -l < gene_count_files.txt)

# if there are at least two files
if [ $n -ge 2 ]; then
    # get Geneid and counts column of the first file
    firstFile=$(head -n 1 gene_count_files.txt)
    tail -n +2 "$firstFile" | cut -f1,7 > file1.txt
    
    # get counts column of remaining files
    i=2
    for file in $(tail -n +2 gene_count_files.txt); do
        tail -n +2 "$file" | cut -f7 > "file$i.txt"
        i=$((i+1))
    done

    # paste files together
    paste -d "\t" file*.txt > gene_counts_matrix.tsv

    # rename columns
    sed -i 's,starAligned/,,g' gene_counts_matrix.tsv
    sed -i 's,.Aligned.sortedByCoord.out.bam,,g' gene_counts_matrix.tsv
    sed -i 's/Geneid/gene_id/g' gene_counts_matrix.tsv

    # cleanup
    rm file*.txt
    
    # print success message
    echo "Gene counts matrix successfully generated."
else
    echo "Could not generate gene counts matrix. You need at least two samples to generate this matrix."
fi

# cleanup
rm gene_count_files.txt


# SUMMED EXONS
#-------------------------------------------------------------------------------

# create list of exon featureCount files
ls -1 | grep exon.counts$ > exon_files.txt

# store number of files
n=$(wc -l < exon_files.txt)

# if there is at least 1 file
if [ $n -ge 1 ]; then
    # get name of one file
    exonFile=$(head -n 1 exon_files.txt)
    
    # store the Geneid and length (summed exons) column
    tail -n +2 "$exonFile" | cut -f1,6 > summed_exons.tsv

    # change Geneid to gene_id and Length to exonic_length
    sed -i 's/Geneid/gene_id/g' summed_exons.tsv
    sed -i 's/Length/exonic_length/g' summed_exons.tsv

    # print success message
    echo "Length of summed exons successfully generated."
else
    echo "Cannot create summed_exons.tsv file. No exon featureCount files exist."
fi

# cleanup
rm exon_files.txt
