#!/bin/bash -l

PINDEL=/apps/CLC_ExternalApps/pindel/pindel
SAMTOOLS=/apps/bio/apps/samtools/1.3.1/samtools
SAM2PINDEL=/apps/CLC_ExternalApps/pindel/sam2pindel
TMPOUT=/medstore/CLC_External_Tmp/pindel_tmp
PINDEL2VCF=/apps/CLC_ExternalApps/pindel/pindel2vcf
DATE=$(date | sed 's/ /_/g' | sed 's/:/_/'g | cut -d"_" -f-6)

while getopts :b:s:f:o:h opt; do
  case $opt in
    b)
        echo "-b (bam) was input as $OPTARG" >&2
        BAM=$OPTARG
    ;;
    s)
        echo "-s (insert size) was input as $OPTARG" >&2
        SIZE=$OPTARG
    ;;
    f)
        echo "-f (reference fasta) was input as $OPTARG" >&2
        FASTA=$OPTARG
    ;;
    o)
        echo "-o (output) was input as $OPTARG" >&2
        OUTPUT=$OPTARG
    ;;
    h)
        echo "$HELP"
        exit 1
    ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        echo "Type $0 -h for usage"
        exit 1
    ;;
  esac
done

BAMFILE=$(basename $BAM)

$SAMTOOLS view $BAM | $SAM2PINDEL - ${TMPOUT}/${BAMFILE} $SIZE sampletag 0 

PINPUT=${TMPOUT}/${BAMFILE}

$PINDEL -f $FASTA -p $PINPUT -c ALL -o ${PINPUT}.pindelout -T 8

$PINDEL2VCF -p ${PINPUT}.pindelout -r $FASTA -R hg19 -d $DATE -v $OUTPUT

wait

# 1rm ${PINPUT}
