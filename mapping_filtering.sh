#!/bin/bash

# Usage: mapping_filter.sh <reference>.fasta <file1>.fastq <file2>.fastq
# Dependencies: bowtie2, samtools, seqtk

echo "Pease input desired name of output directory: "
read OUTDIR

echo `hostname`
date

DBNAME=`date +%C%y_%m_%d`
SAM_FULL=mapping_full_$DBNAME.sam
SAM_MAPPER=mappers_$DBNAME.sam
SAM_NON_MAPPER=non_mappers_$DBNAME.sam
LIST_MAPPER=mappers_$DBNAME.lst
LIST_NON_MAPPER=non_mappers_$DBNAME.lst

echo "Creating bowtie2 database"

# Creates a bowtie2 database and names it by date and a random number
if [ ! -e $DBNAME ]; then

	bowtie2-build -f $1 $DBNAME
else
	echo "DB exists, proceeding..." 
fi

wait

mkdir $OUTDIR

# Starts bowtie2 mapping
bowtie2 -x $DBNAME -1 $2 -2 $3 -S $OUTDIR/$SAM_FULL 

# Splits the sam files into mappers and non_mappers
samtools view -S -F4 $OUTDIR/$SAM_FULL > $OUTDIR/$SAM_MAPPER &
samtools view -S -f4 $OUTDIR/$SAM_FULL > $OUTDIR/$SAM_NON_MAPPER &
wait

#Makes lists containing the headers of the mapping and non_mapping reads
cut -f1 $OUTDIR/$SAM_MAPPER | sort | uniq > $OUTDIR/$LIST_MAPPER &
cut -f1 $OUTDIR/$SAM_NON_MAPPER | sort | uniq > $OUTDIR/$LIST_NON_MAPPER &

seqtk subseq $1 $OUTDIR/$LIST_MAPPER > $OUTDIR/mappers_$2 &
seqtk subseq $2 $OUTDIR/$LIST_MAPPER > $OUTDIR/mappers_$3 &
wait

seqtk subseq $1 $OUTDIR/$LIST_NON_MAPPER > $OUTDIR/non_mappers_$2 &
seqtk subseq $2 $OUTDIR/$LIST_NON_MAPPER > $OUTDIR/non_mappers_$3 &
wait

date
echo "Finished running bowtie2 mapping"
