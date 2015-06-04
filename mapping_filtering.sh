#!/bin/bash

# Usage: mapping_filter.sh <reference>.fasta <file1>.fastq <file2>.fastq
# Dependencies: bowtie2, samtools, seqtk

echo "Pease input desired name of your analysis: "
read OUTDIR

echo "Do you want the intermediary files to be deleted after completion?(y/n)"
read DELETE

echo `hostname`
date

DATE=`date +%C%y_%m_%d`
SAM_FULL=mapping_full_$DATE.sam
SAM_MAPPER=mappers_$DATE.sam
SAM_NON_MAPPER=non_mappers_$DATE.sam
LIST_MAPPER=mappers_$DATE.lst
LIST_NON_MAPPER=non_mappers_$DATE.lst

echo "Creating bowtie2 database..."
# Creates a bowtie2 database and names it by date and a random number
if [ ! -e $1.?.bt2 ]; then

	bowtie2-build -f $1 $1
else
	echo "Database aready exists, proceeding..." 
fi

wait

mkdir $OUTDIR

# Starts bowtie2 mapping
echo "Running bowtie2 mapping..."
bowtie2 -x $1 -1 $2 -2 $3 -S $OUTDIR/$SAM_FULL 
wait

# Splits the sam files into mappers and non_mappers
echo "Separating sam files..."
samtools view -S -F4 $OUTDIR/$SAM_FULL > $OUTDIR/$SAM_MAPPER &
samtools view -S -f4 $OUTDIR/$SAM_FULL > $OUTDIR/$SAM_NON_MAPPER &
wait

mkdir $OUTDIR/lists
mkdir $OUTDIR/mapped_reads
mkdir $OUTDIR/non_mapped_reads

# Makes lists containing the headers of the mapping and non_mapping reads
echo "Creating lists..."
cut -f1 $OUTDIR/$SAM_MAPPER | sort | uniq > $OUTDIR/lists/"$OUTDIR"_$LIST_MAPPER &
cut -f1 $OUTDIR/$SAM_NON_MAPPER | sort | uniq > $OUTDIR/lists/"$OUTDIR"_$LIST_NON_MAPPER &
wait

# Pulling mapped reads from libraries
echo "Fetching reads..."
seqtk subseq $2 $OUTDIR/lists/"$OUTDIR"_$LIST_MAPPER > $OUTDIR/mapped_reads/"$OUTDIR"_"$DATE"_mappers_$2 &
seqtk subseq $3 $OUTDIR/lists/"$OUTDIR"_$LIST_MAPPER > $OUTDIR/mapped_reads/"$OUTDIR"_"$DATE"_mappers_$3 &
wait

# Pulling non-mapped reads from libraries
seqtk subseq $2 $OUTDIR/lists/"$OUTDIR"_$LIST_NON_MAPPER > $OUTDIR/non_mapped_reads/"$OUTDIR"_"$DATE"_non_mappers_$2 &
seqtk subseq $3 $OUTDIR/lists/"$OUTDIR"_$LIST_NON_MAPPER > $OUTDIR/non_mapped_reads/"$OUTDIR"_"$DATE"_non_mappers_$3 &
wait

# Deletes intermediary files if desired
if [[ $DELETE = y* ]] || [[ $DELETE = Y* ]]; then
	echo "Removing intermediary sam files..."
	rm $OUTDIR/*.sam
else
	echo "Keeping sam files"
fi

date
echo "Finished running bowtie2 mapping and filtering"
