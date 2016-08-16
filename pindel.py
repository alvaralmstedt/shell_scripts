#!/bin/bash -l

PINDEL=/apps/CLC_ExternalApps/pindel/pindel
SAMTOOLS=/apps/bio/apps/samtools/1.3.1/samtools
SAM2PINDEL=/apps/CLC_ExternalApps/pindel/sam2pindel
TMPOUT=/apps/CLC_Tmp

while getopts :d:s:f:i:n:kht: opt; do
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
    n)
        echo "-n (name) was input as $OPTARG" >&2
        NAME=$OPTARG
    ;;
    k)
        echo "-k (keep) was triggered, sam files will be kept" >&2
        DELETE=0
    ;;
    h)
        echo "$HELP"
        exit 1
    ;;
    t)
        echo "-t number of cores was input as $OPTARG" >&2
        CPU=$OPTARG
    ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        echo "Type $0 -h for usage"
        exit 1
    ;;
  esac
done

samtools view $BAM | $SAM2PINDEL - ${TMPOUT}/output4pindel.txt $SIZE sampletag 0 

$PINDEL -f $FASTA -i $CONFIG -c ALL -o ${TMPOUT}/${BAM}.pindelout -T 8
