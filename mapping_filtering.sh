#!/bin/bash

# Usage: mapping_filter.sh <reference>.fasta <file1>.fastq <file2>.fastq
# Dependencies: bowtie2, samtools, seqtk

echo "Name of outputdirectory: "
read OUTDIR

echo `hostname`
date

BNAME=`date +%C%y_%m_%d`
SAM_FULL=mapping_full_$DBNAME.sam
SAM_MAPPER=mappers_$DBNAME.sam
SAM_NON_MAPPER=non_mappers_$DBNAME.sam
LIST_MAPPER=mappers_$DBNAME.lst
LIST_NON_MAPPER=non_mappers_$DBNAME.lst

echo "Creating bowtie2 database"

# Creates a bowtie2 database and names it by date and a random number
if [ ! -e $DBNAME ];
do
	bowtie2-build -f $1 $DBNAME 
done

# Starts bowtie2 mapping
bowtie2 -x $DBNAME -1 $1 -2 $2 -S $OUTDIR/$SAM_FULL 

# Splits the sam files into mappers and non_mappers
samtools view -S F4 $OUTDIR/$SAM_FULL > $OUTDIR/$SAM_MAPPER &
samtools view -S f4 $OUTDIR/$SAM_FULL > $OUTDIR/$SAM_NON_MAPEPR &
wait

#Makes lists containing the headers of the mapping and non_mapping reads
cut -f1 $OUTDIR/$SAM_MAPPER | sort | uniq > $OUTDIR/$LIST_MAPPER &
cut -f1 $OUTDIR/$SAM_NON_MAPPER | sort | uniq > $OUTDIR/$LIST_NON_MAPPER &

seqtk subseq $1 $OUTIDR/$MAPPER_LIST > $OUTIDR/mappers_$1 &
seqtk subseq $2 $OUTIDR/$MAPPER_LIST > $OUTIDR/mappers_$2 &
wait

seqtk subseq $1 $OUTIDR/$NON_MAPPER_LIST > $OUTIDR/non_mappers_$1 &
seqtk subseq $2 $OUTIDR/$NON_MAPPER_LIST > $OUTIDR/non_mappers_$2 &
wait

date
echo "Finished running bowtie2 mapping"
