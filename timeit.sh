#!/bin/bash -l

INFILE=$1

module load samtools/1.3.1

COMMAND_PIPE="samtools sort $INFILE -n -@ 40 -m 2G -T /tmp/timetest | samtools view - -o /tmp/timetest/output_pipe.sam -@ 40 -h"

COMMAND_AMP="samtools sort $INFILE -n -@ 40 -m 2G -T /tmp/timetest -o /tmp/timetest/output.bam && samtools view /tmp/timetest/output.bam -@ 40 -h -o /tmp/timetest/output_amp.sam"

{ time (eval $COMMAND_PIPE) ; } 2> /tmp/timetest/pipetime.txt
wait
ls -l /tmp/timetest/output_pipe.sam >> /tmp/timetest/pipetime.txt
rm /tmp/timetest/output_pipe.sam
wait
{ time (eval $COMMAND_AMP) ; } 2> /tmp/timetest/amptime.txt
wait
ls -l /tmp/timetest/output_amp.sam >> /tmp/timetest/amptime.txt
rm /tmp/timetest/output.bam
rm /tmp/timetest/output_amp.sam
