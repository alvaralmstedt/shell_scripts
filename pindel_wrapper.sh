#!/bin/bash -l

set -e

PINDEL=/apps/CLC_ExternalApps/pindel/pindel
SAMTOOLS=/apps/bio/apps/samtools/1.3.1/samtools
SAM2PINDEL=/apps/CLC_ExternalApps/pindel/sam2pindel
TMPOUT=/medstore/CLC_External_Tmp/pindel_tmp
PINDEL2VCF=/apps/CLC_ExternalApps/pindel/pindel2vcf
DATE=$(date | sed 's/ /_/g' | sed 's/:/_/'g | cut -d"_" -f-6)

while getopts :b:s:f:o:h opt; do
  case $opt in
    b)
        echo "-b (bam) was input as $OPTARG" >&1
        BAM=$OPTARG
    ;;
    s)
        echo "-s (insert size) was input as $OPTARG" >&1
        SIZE=$OPTARG
    ;;
    f)
        echo "-f (reference fasta) was input as $OPTARG" >&1
        FASTA=$OPTARG
    ;;
    o)
        echo "-o (output) was input as $OPTARG" >&1
        OUTPUT=$OPTARG
    ;;
    h)
        echo "$HELP"
        exit 1
    ;;
    \?)
        echo "Invalid option: -$OPTARG" >&1
        echo "Type $0 -h for usage"
        exit 1
    ;;
  esac
done

BAMFILE=$(basename $BAM)

if [[ $FASTA == *"hg19"* ]] ; then
    echo "Using default hg19"
    FASTA=/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta
else
    echo "Running samtools faidx on $FASTA"
    $SAMTOOLS faidx $FASTA
fi

echo "Samtools faidx completed"
echo "Starting samtools view + sam2pindel"

$SAMTOOLS view $BAM | $SAM2PINDEL - ${TMPOUT}/${BAMFILE}.txt $SIZE sampletag 0 Illumina-PairEnd

echo "Samtools view + sam2pindel completed"

PINPUT=${TMPOUT}/${BAMFILE}.txt

echo "Starting pindel"

$PINDEL -f $FASTA -p $PINPUT -c ALL -o ${PINPUT}.pindelout -T 8

echo "Pindel completed, concatenationg results"

for i in $(ls ${TMPOUT} | grep -v "concatenate"); do
    cat ${TMPOUT}/$i >> ${TMPOUT}/concatenatedpindel.out
done

echo "Results condatenated, starting pindel2vcf"

$PINDEL2VCF -p ${TMPOUT}/concatenatedpindel.out -r $FASTA -R hg19 -d $DATE -v $OUTPUT

echo "Finished"
