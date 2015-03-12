#!/bin/bash

# Note: contig_average_coverage.py and fasta_analyzer.py have to be added to your $PATH for this
# script to work. Courtesy to Sandra Karlsten for these python programs.


# 2 arguments: smrtanalysis job number (ex 016536) and output prefix/name (whatever you want)
SMRT_JOB_NUM=$1
OUT_PREFIX=$2

echo "Creating symlink and copying assembly..."
# Create a symlink to the bam-file in the selected job numer folder on DNA server.
ln -s /data01/smrtanalysis/userdata/jobs/016/$SMRT_JOB_NUM/data/aligned_reads.bam
cp /data01/smrtanalysis/userdata/jobs/016/$SMRT_JOB_NUM/data/polished_assembly.fasta.gz ./
echo "Symlink created and assembly copied..."

gunzip polished_assembly.fasta.gz
echo "Assembly gunzipped..."

echo "Commencing generation of coverage file..."
# Generates a coverage file from
genomeCoverageBed -ibam aligned_reads.bam -d > $OUT_PREFIX.coverage
echo "Coverage file generated..."

echo "Cutting..."
# Cuts out only the contig name and cov per base value, which are what we want
cut -f 1,3 $OUT_PREFIX.coverage > cut_$OUT_PREFIX.coverage
echo "Cut made..."
echo "Launching contig_average_coverage.py..."
#Calculates and replaces per base coverage with avg cov for the contig
contig_average_coverage.py cut_$OUT_PREFIX.coverage > cut_$OUT_PREFIX.avg.coverage
echo "Average calculated..."

echo "Deleting..."
# Remove intermediary files
rm $OUT_PREFIX.coverage
rm cut_$OUT_PREFIX.coverage
echo "Intermediary files removed..."

echo "Formatting fasta file headers..."
sed 's/|quiver//g' polished_assembly.fasta > temp && mv temp polished_assembly.fasta
echo "Headers formatted..."

echo "Launching fasta_analyzer.py"
# Launch fasta_analyzer.py in coverage/gc plot mode
fasta_analyzer.py -cg polished_assembly.fasta cut_$OUT_PREFIX.avg.coverage
