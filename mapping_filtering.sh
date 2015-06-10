#!/bin/bash

# Dependencies: bowtie2, samtools, seqtk
# Some code borrowed from https://wikis.utexas.edu/display/bioiteam/Example+BWA+alignment+script

# Test if the correct number of input files are specified

if [ "$2" == "" ]; then

	echo ""
	echo "Usage: mapping_filter.sh <reference>.fasta <file1>.fastq [<file2>.fastq]"
	echo ""
	echo "	A reference fasta file and one (for singlets)" 
	echo "	or two (paired end) fastq files are required."
	echo ""
	exit 1;
fi

	echo -n "Name of output directory: "
	read NAME

	echo -n "Remove *.sam files after completion? (y/n): "
	read DELETE

# Print some informative error meassages
err() {
    echo "$1 exited unexpectedly";
    exit 1;
}

# Function for checking the exit code of a child process

ckeckExit() {
    if [ "$1" == "0" ]; then
        echo "[Done] $2 `date`";
    else
        err "[Error] $2 returned non-0 exit code $1";
    fi
}

DATE=`date +%C%y_%m_%d`
SAM_FULL=mapping_full_$DATE.sam
SAM_MAPPER=mappers_$DATE.sam
SAM_NON_MAPPER=non_mappers_$DATE.sam
LIST_MAPPER=mappers_$DATE.lst
LIST_NON_MAPPER=non_mappers_$DATE.lst
LIST_TRUE_NON_MAPPER=non_mapper.lst
LIST_TRUE_MAPPER=mapper.lst
OUTDIR="$NAME"_$DATE
MAPPING_INFO=README_mapping.txt

# Creates a bowtie2 database and names it by date and a random number

files=$(ls "$1".?.bt2 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then
	echo "[info] Creating bowtie2 database..."
	bowtie2-build -f $1 $1
	ckeckExit $? "bowtie2-build"
else
	echo "[info] Database aready exists, proceeding..." 
fi

wait

mkdir $OUTDIR

# Starts bowtie2 mapping
echo "[info] Running bowtie2 mapping..."
bowtie2 -x $1 -1 $2 -2 $3 -S $OUTDIR/$SAM_FULL 2> $OUTDIR/$MAPPING_INFO
ckeckExit $? "bowtie2"
wait

# Splits the sam files into mappers and non_mappers
echo "[info] Separating sam files..."
samtools view -S -F4 $OUTDIR/$SAM_FULL > $OUTDIR/$SAM_MAPPER &
    ckeckExit $? "samtools"
samtools view -S -f4 $OUTDIR/$SAM_FULL > $OUTDIR/$SAM_NON_MAPPER &
    ckeckExit $? "samtools"
wait

mkdir $OUTDIR/lists
mkdir $OUTDIR/mapped_reads
mkdir $OUTDIR/non_mapped_reads
mkdir $OUTDIR/half_mapped_reads

# Makes lists containing the headers of the mapping and non_mapping reads
echo "[info] Creating lists..."
cut -f1 $OUTDIR/$SAM_MAPPER | sort | uniq > $OUTDIR/lists/"$NAME"_$LIST_MAPPER &
    ckeckExit $? "cut"
cut -f1 $OUTDIR/$SAM_NON_MAPPER | sort | uniq > $OUTDIR/lists/"$NAME"_$LIST_NON_MAPPER &
    ckeckExit $? "cut"
wait

diff $OUTDIR/lists/"$NAME"_$LIST_NON_MAPPER $OUTDIR/lists/"$NAME"_$LIST_MAPPER | grep "> " | sed "s/> //g" > $OUTDIR/lists/"$NAME"_$LIST_TRUE_MAPPER
    ckeckExit $? "diff, grep or sed"
wait

diff $OUTDIR/lists/"$NAME"_$LIST_NON_MAPPER $OUTDIR/lists/"$NAME"_$LIST_MAPPER | grep "< " | sed "s/< //g" > $OUTDIR/lists/"$NAME"_$LIST_TRUE_NON_MAPPER
    ckeckExit $? "diff, grep or sed"
wait

touch $OUTDIR/lists/temp.lst

cat $OUTDIR/lists/"$NAME"_$LIST_MAPPER > $OUTDIR/lists/temp.lst
    ckeckExit $? "cat"
cat $OUTDIR/lists/"$NAME"_$LIST_NON_MAPPER >> $OUTDIR/lists/temp.lst
    ckeckExit $? "cat"
wait

sort $OUTDIR/lists/temp.lst | uniq -d > $OUTDIR/lists/"$NAME"_half_mappers.lst
    ckeckExit $? "sort"
wait

# Comment the following three lines out if you want to double-check list numbers
rm $OUTDIR/lists/temp.lst
rm $OUTDIR/lists/"$NAME"_$LIST_MAPPER
rm $OUTDIR/lists/"$NAME"_$LIST_NON_MAPPER
wait

# Pulling mapped reads from libraries
echo "[info] Fetching reads..."
# seqtk subseq $2 $OUTDIR/lists/"$NAME"_$LIST_MAPPER > $OUTDIR/mapped_reads/"$NAME"_mappers_$2 &
# seqtk subseq $3 $OUTDIR/lists/"$NAME"_$LIST_MAPPER > $OUTDIR/mapped_reads/"$NAME"_mappers_$3 &
# wait

seqtk subseq $2 $OUTDIR/lists/"$NAME"_$LIST_TRUE_MAPPER > $OUTDIR/mapped_reads/"$NAME"mappers_$2 &
    ckeckExit $? "seqtk"
seqtk subseq $3 $OUTDIR/lists/"$NAME"_$LIST_TRUE_MAPPER > $OUTDIR/mapped_reads/"$NAME"mappers_$3 &
    ckeckExit $? "seqtk"
wait

# Pulling non-mapped reads from libraries
# seqtk subseq $2 $OUTDIR/lists/"$NAME"_$LIST_NON_MAPPER > $OUTDIR/non_mapped_reads/"$NAME"_non_mappers_$2 &
# seqtk subseq $3 $OUTDIR/lists/"$NAME"_$LIST_NON_MAPPER > $OUTDIR/non_mapped_reads/"$NAME"_non_mappers_$3 &
# wait

seqtk subseq $2 $OUTDIR/lists/"$NAME"_$LIST_TRUE_NON_MAPPER > $OUTDIR/non_mapped_reads/"$NAME"non_mappers_$2 &
    ckeckExit $? "seqtk"
seqtk subseq $3 $OUTDIR/lists/"$NAME"_$LIST_TRUE_NON_MAPPER > $OUTDIR/non_mapped_reads/"$NAME"non_mappers_$3 &
    ckeckExit $? "seqtk"
wait

# Pulling half_mapped reads from libraries
seqtk subseq $2 $OUTDIR/lists/"$NAME"_half_mappers.lst > $OUTDIR/half_mapped_reads/"$NAME"half_mappers_$2 &
    ckeckExit $? "seqtk"
seqtk subseq $3 $OUTDIR/lists/"$NAME"_half_mappers.lst > $OUTDIR/half_mapped_reads/"$NAME"half_mappers_$3 &
    ckeckExit $? "seqtk"
wait

# Deletes intermediary files if desired
if [[ $DELETE = y* ]] || [[ $DELETE = Y* ]]; then
	echo "[info] Removing intermediary sam files..."
	rm $OUTDIR/*.sam
    else
	echo "[info] Keeping sam files"
fi

echo "Finished running bowtie2 mapping and filtering $(date)"
