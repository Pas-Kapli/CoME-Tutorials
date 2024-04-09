#!/bin/bash

# This is so you can run the script multiple times
# with adding the sequences multiple times
rm -f msa.txt

for seqfile in $(cat sequenceFiles.txt)
do
	name=$(cat $seqfile | awk '{ print $1}')

	# Check to make sure sequence isn't empty
	if [[  $(cat $seqfile | awk '{print $2}' | grep [^-]) ]]; then

		seq=$(cat $seqfile | awk '{print $2}')

	fi

	# Add sequence to file
	echo ${name} >> msa.txt
	echo ${seq} >> msa.txt


done


raxml-ng --msa msa.txt --model JC --msa-format FASTA
