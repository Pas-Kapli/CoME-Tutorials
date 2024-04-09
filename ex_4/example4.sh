#!/bin/bash

# This is so the script can be run multiple times
# without needing to remove the sequence file 
# created previously.
rm -f sequenceToAlign.txt

for seq in $(cat sequenceFiles.txt)
do

	# Add sequence to file
	cat ${seq} >> sequenceToAlign.txt

done

# Formatting the sequences for muscle
sed -i 's/\t/\n/g' sequenceToAlign.txt

muscle -align sequenceToAlign.txt -output msa_align.fa &> outMuscle 

