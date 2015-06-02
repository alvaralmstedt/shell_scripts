#!/bin/bash
echo `hostname`
date

#DBNAME='mapping_18052015'
#READPATH='/data01/alvar/mapping_11th_may_2015/pulled/non_mappers'
#MPPATH="/data02/tomas/skeletonema_MP"
#PEPATH="/data01/alvar/mapping_4th_may_2015/150_300_650_pair"

echo "Running bowtie2 mapping"

bowtie2 -x mapping_01062015 -1 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/P1354_Fw_quad_non_mappers.fastq -2 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/P1354_Rev_quad_non_mappers.fastq -S results/P1354_mapping_01062015.sam &

bowtie2 -x mapping_01062015 -1 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/P1355_Fw_quad_non_mappers.fastq -2 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/P1355_Rev_quad_non_mappers.fastq -S results/P1355_mapping_01062015.sam &

bowtie2 -x mapping_01062015 -1 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/150Pair_Fw_quad_non_mappers.fastq -2 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/150Pair_Rev_quad_non_mappers.fastq -S results/150Pair_mapping_01062015.sam &

bowtie2 -x mapping_01062015 -1 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/300Pair_Fw_quad_non_mappers.fastq -2 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/300Pair_Rev_quad_non_mappers.fastq -S results/300Pair_mapping_01062015.sam &

bowtie2 -x mapping_01062015 -1 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/650Pair_Fw_quad_non_mappers.fastq -2 /home/alvar/mapping_25th_may_2015/pulled/non_mappers/650Pair_Rev_quad_non_mappers.fastq -S results/650Pair_mapping_01062015.sam &


date
echo "Finished running bowtie2 mapping"
