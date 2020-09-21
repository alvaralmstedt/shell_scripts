#!/bin/bash -l

TMPDIR=/tmp/svbay
DATADIR=$TMPDIR/sv-bay-data
BAMDIR=$DATADIR/bam
FASTADIR=$DATADIR/fa_files
GEMDIR=$DATADIR/gem_files
SAMTOOLS=/apps/bio/apps/samtools/1.3.1/samtools
SEPSAM=/apps/CLC_ExternalApps/SV-Bay/SV-Bay/utils/separately_save_sam2.py
FREEC=/apps/CLC_ExternalApps/SV-Bay/SV-Bay/Control-FREEC/freec
FREECCONFIG1=/apps/CLC_ExternalApps/SV-Bay/SV-Bay/Control-FREEC/configs/config1.txt
FREECCONFIG2=/apps/CLC_ExternalApps/SV-Bay/SV-Bay/Control-FREEC/configs/config2.txt
BAY_TUM_CONFIG=/apps/CLC_ExternalApps/SV-Bay/SV-Bay/config/config.yaml
BAY_NORM_CONFIG=/apps/CLC_ExternalApps/SV-Bay/SV-Bay/config/config_germ.yaml
TEXTINPUT=False

while getopts :b:c:n:m:o:h opt; do
  case $opt in
    b)
        echo "-b (bam (tumour)) was input as $OPTARG" >&1
        BAM=$OPTARG
    ;;
    c)
        echo "-c (tumour text) was input as $OPTARG" >&1
        TUMOUR_TEXT=$OPTARG
    ;;    
    n)
        echo "-n (normal bam) was input as $OPTARG" >&1
        NORMAL=$OPTARG
    ;;
    m)
        echo "-m (normal text) was input as $OPTARG" >&1
        NORMAL_TEXT=$OPTARG
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

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load anaconda2/4.1.0
module load samtools/1.3.1

###################### SETUP #######################

if [ ! -z $TUMOUR_TEXT  ] && [ ! -z $NORMAL_TEXT ] ; then
    BAM=$(cat $TUMOUR_TEXT)
    NORMAL=$(cat $NORMAL_TEXT)
    TEXTINPUT=True
elif [ ! -z $TUMOUR_TEXT ] ; then
    BAM=$(cat $TUMOUR_TEXT)
elif [ ! -z $NORMAL_TEXT ] ; then
    NORMAL=$(cat $NORMAL_TEXT)
fi

if [[ ! -z $NORMAL ]] ; then

DATADIR2=$TMPDIR/sv-bay-data2
BAMDIR2=$DATADIR2/bam
FASTADIR2=$DATADIR2/fa_files
GEMDIR2=$DATADIR2/gem_files

echo "HOSTNAME = $HOSTNAME"
echo "BAM = $BAM"
echo "NORMAL = $NORMAL"
echo "OUTPUT = $OUTPUT"

mkdir $TMPDIR
mkdir $TMPDIR/tumour
mkdir $TMPDIR/normal
mkdir $DATADIR
mkdir $DATADIR2
mkdir $BAMDIR
mkdir $BAMDIR2
mkdir $FASTADIR
mkdir $FASTADIR2
mkdir $GEMDIR
mkdir $GEMDIR2

echo "Directories created at `date`"

cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/chromosomes/* $FASTADIR/
cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/chromosomes/* $FASTADIR2/
cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/gem_hg19/* $GEMDIR/
cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/gem_hg19/* $GEMDIR2/
cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/centrom_hg19.txt $DATADIR/
cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/centrom_hg19.txt $DATADIR2/

else

echo "HOSTNAME = $HOSTNAME"
echo "BAM = $BAM"
echo "OUTPUT = $OUTPUT"

mkdir $TMPDIR
mkdir $DATADIR
mkdir $BAMDIR
mkdir $FASTADIR
mkdir $GEMDIR

echo "Directories created at `date`"

cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/chromosomes/* $FASTADIR/
cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/gem_hg19/* $GEMDIR/
cp /apps/CLC_ExternalApps/SV-Bay/SV-Bay/data/centrom_hg19.txt $DATADIR/
fi

echo "Files copied at `date`"

#${SAMTOOLS} sort $BAM -o ${TMPDIR}/chr$(basename ${BAM%.*})_sorted.bam -@ 5 -m 2G >&1

#################### FREEC ######################

#if [[ ! -z $NORMAL]] ; then
#    $FREEC -conf $FREECCONFIG1
#    $FREEC -conf $FREECCONFIG2
#else
#    $FREEC -conf $FREECCONFIG1
#fi
##################### SVBAY #####################

if [[ ! -z $NORMAL ]] ; then
    #------------------------------------------------------
    #------------------TUMOUR NORMAL-----------------------
    #------------------------------------------------------

#    if [ $TEXTINPUT = False ]; then
#    ${SAMTOOLS} sort $BAM -o ${TMPDIR}/tumour/chr$(basename ${BAM%.*})_sorted.bam -@ 5 -m 2G & >&1
cp $BAM $TMPDIR/tumour/tumour.bam
${SAMTOOLS} sort $TMPDIR/tumour/tumour.bam -o ${TMPDIR}/tumour/tumour_sorted.bam -@ 40 -m 1G >&1
#    ${SAMTOOLS} sort $NORMAL -o ${TMPDIR}/normal/chr$(basename ${NORMAL%.*})_sorted.bam -@ 5 -m 2G >&1
echo "Tumour bam sorted at `date`"

cp $NORMAL $TMPDIR/normal/normal.bam
${SAMTOOLS} sort $TMPDIR/normal/normal.bam -o ${TMPDIR}/normal/normal_sorted.bam -@ 20 -m 2G >&1
echo "Normal bam sorted at `date`"
#    else
#        cp $BAM $TMPDIR/tumour/tumour.bam
#        cp $NORMAL $TMPDIR/normal/normal.bam
#        echo "Bam files copied over at `date`"
#    fi     

    $FREEC -conf $FREECCONFIG1 >&1
    $FREEC -conf $FREECCONFIG2 >&1
    echo "Freec runs completed at `date`"

#    python $SEPSAM -i ${TMPDIR}/tumour/chr$(basename ${BAM%.*})_sorted.bam -o $TMPDIR/tumour/ >&1
    python $SEPSAM -i ${TMPDIR}/tumour/tumour_sorted.bam -o $TMPDIR/tumour/ >&1
#    python $SEPSAM -i ${TMPDIR}/normal/chr$(basename ${NORMAL%.*})_sorted.bam -o $TMPDIR/normal/ >&1
    python $SEPSAM -i ${TMPDIR}/normal/normal_sorted.bam -o $TMPDIR/normal/ >&1

    echo "Bams separated at `date`"

    for i in $TMPDIR/tumour/*.sam ; do
        cat $i | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0}  $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > ${BAMDIR}/chr$(basename ${i%.*})_sorted.bam >&1
    done

    for i in $TMPDIR/normal/*.sam ; do
        cat $i | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0}  $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > ${BAMDIR2}/chr$(basename ${i%.*})_sorted.bam >&1
    done

    echo "Sams renamed and converted into bam at `date`"

    echo "ls TMPDIR:"
    ls $TMPDIR >&1

    for i in ${BAMDIR}/*.?am ; do
        ${SAMTOOLS} index ${BAMDIR}/$(basename ${i%.*}).bam >&1
    done

    for i in ${BAMDIR2}/*.?am ; do
        ${SAMTOOLS} index ${BAMDIR2}/$(basename ${i%.*}).bam >&1
    done

    echo "Bams indexed at `date`"

    echo "ls BAMDIR($BAMDIR):"
    ls $BAMDIR >&1
    echo "ls BAMDIR2($BAMDIR2):"
    ls $BAMDIR2 >&1

    echo "Bams sorted and indexed at `date`"

    python -B /apps/CLC_ExternalApps/SV-Bay/SV-Bay/src/main_clustering.py -c $BAY_TUM_CONFIG &
    python -B /apps/CLC_ExternalApps/SV-Bay/SV-Bay/src/main_clustering.py -c $BAY_NORM_CONFIG &
    python -B /apps/CLC_ExternalApps/SV-Bay/SV-Bay/src/main_probabilities.py -c $BAY_TUM_CONFIG
    python -B /apps/CLC_ExternalApps/SV-Bay/SV-Bay/src/main_assemly_links.py -c $BAY_TUM_CONFIG -n "$DATADIR2/cluster_files" > $DATADIR/svbay_results.txt

    echo "SV-bay runs completed at `date`"

else    
    #------------------------------------------------
    #-----------------SINGLE SAMPLE------------------
    #------------------------------------------------
    

#    if [ $TEXTINPUT = False ]; then
        cp $BAM $TMPDIR/tumour.bam
        ${SAMTOOLS} sort $TMPDIR/tumour.bam -o ${TMPDIR}/tumour_sorted.bam -@ 40 -m 1G >&1
#        echo "Bams sorted at `date`"
#    else
#    fi
    
    echo "Bams sorted at `date`"
    
    $FREEC -conf $FREECCONFIG1 >&1
    echo "Freec completed at `date`"

    python $SEPSAM -i ${TMPDIR}/tumour_sorted.bam -o $TMPDIR/ >&1
#    python $SEPSAM -i ${TMPDIR}/chr$(basename ${BAM%.*})_sorted.bam -o $TMPDIR/ >&1

    echo "Bams separated at `date`"

    for i in $TMPDIR/*.sam ; do
        cat $i | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0}  $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > ${BAMDIR}/chr$(basename ${i%.*})_sorted.bam >&1
    done

    echo "Sams renamed and converted into bam at `date`"

    echo "ls TMPDIR"
    ls $TMPDIR >&1

    for i in ${BAMDIR}/*.?am ; do
#    ${SAMTOOLS} sort $i -o ${BAMDIR}/chr$(basename ${i%.*})_sorted.bam >&1
        ${SAMTOOLS} index ${BAMDIR}/$(basename ${i%.*}).bam >&1
    done

    echo "Bams indexed at `date`"

    echo "ls BAMDIR"
    ls $BAMDIR >&1
    echo "Bams sorted and indexed at `date`"

    python -B /apps/CLC_ExternalApps/SV-Bay/SV-Bay/src/main_clustering.py -c /apps/CLC_ExternalApps/SV-Bay/SV-Bay/config/config.yaml
    python -B /apps/CLC_ExternalApps/SV-Bay/SV-Bay/src/main_probabilities.py -c /apps/CLC_ExternalApps/SV-Bay/SV-Bay/config/config.yaml
    python -B /apps/CLC_ExternalApps/SV-Bay/SV-Bay/src/main_assemly_links.py -c /apps/CLC_ExternalApps/SV-Bay/SV-Bay/config/config.yaml > $DATADIR/svbay_results.txt

    echo "SV-bay runs completed at `date`"

fi

cp $DATADIR/svbay_results.txt $OUTPUT

echo "Data moved at `date`"

#rm -rf $TMPDIR

#echo "temp-files deleted"
