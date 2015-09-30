#!/bin/bash

# This script is for concatenating multiple aligned fasta files into one nexus file
# Dependencies: names_first.py, fasta2tab, tab2fasta, seqret

# Run from cwd

while getopts :s: opt; do
  case $opt in
	s)
		echo "-s (your Species) was input as $OPTARG" >&2
		SPECIES=$OPTARG
	;;
  esac
done

if [ $SPECIES = 0 ] ; then
	echo "You must input your species name as an argument"
	exit 1
fi

CWD=$(pwd -P)
FILES=$(ls)
wait
mkdir named
for f in $FILES
do
	echo "Reformatting header of file: $f" 
	names_first.py $f > named/named_$f
done
wait
NAMEDFILES=$(cd named && ls)
wait
mkdir tabbed
echo $NAMEDFILES

for nf in $NAMEDFILES
do
	echo "Tabbing file: $nf"
	fasta2tab $CWD/named/$nf 2> /dev/null | sort > $CWD/tabbed/tabbed_$nf
done
wait

TABBEDFILES=$(cd tabbed && ls)

for spec in $TABBEDFILES
do
	echo $SPECIES
	echo -n "${SPECIES}|gi|" | cat - $CWD/tabbed/$spec | sed '/^$/d' > $CWD/tmp && mv $CWD/tmp $CWD/tabbed/${SPECIES}_${spec}
	head $CWD/tabbed/${SPECIES}_${spec}
	rm $CWD/tabbed/$spec
#	echo "end of specloop"
#	sleep 1
#	sed -i '1i/^/$SPECIES gi|;' spec
done


set -- $TABBEDFILES
wait
mkdir completed
echo $1
sleep 5
# This loop created the files with the correct species names
for header in $CWD/tabbed/$1
do
	echo header
	FILNAMN=$( ( head $header | cut -f1 -d'gi|' ) )
	echo $FILNAMN
	touch $CWD/completed/completed_$FILNAMN
#	cut -f1 -d'     ' header > $CWD/completed/completed_$header
	sleep 1
done
sleep 5
FINALFILES=$(cd completed && ls)
set -- $FINALFILES
INCREMENTER=1

# This loop pastes the species headers into the files
for tf in $TABBEDFILES
do
	echo $1
	cut -f1 -d'	' tf >> $CWD/completed/completed_$(${INCREMENTER})
#	This loop pastes the sequences into the result files	
	for seq in $TABBEDFILES
	do
		cut -f2 -d'	' seq >> $CWD/completed/completed_$($INCREMENTER)
	done
	INCREMENTER=$((INCREMENTER+1))
	sleep 1
done
sleep 5
echo "done"	
