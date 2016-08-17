#!/bin/bash -l

PINDEL=/apps/CLC_ExternalApps/pindel/pindel
SAMTOOLS=/apps/bio/apps/samtools/1.3.1/samtools
SAM2PINDEL=/apps/CLC_ExternalApps/pindel/sam2pindel
TMPOUT=/medstore/CLC_External_Tmp/pindel_tmp

while getopts :d:s:f:i:h opt; do
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
    i)
        echo "-i (config) was input as $OPTARG" >&2
        CONFIG=$OPTARG
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

$SAMTOOLS view $BAM | $SAM2PINDEL - ${TMPOUT}/output4pindel.txt $SIZE sampletag 0 

$PINDEL -f $FASTA -i $CONFIG -c ALL -o ${TMPOUT}/${BAM}.pindelout -T 8

rm ${TMPOUT}/output4pindel.txt
