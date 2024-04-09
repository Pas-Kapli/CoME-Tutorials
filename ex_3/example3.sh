#!/bin/bash

# This is included so you don't duplicate all the sequences
# if you run this script multiple times. 

rm -f sequenceToAlign.txt 

for seq in $(cat sequenceFiles.txt)
do
	# Add sequence to file
	cat ${seq} >> sequenceToAlign.txt

done

sed -i 's/\t/\n/g' sequenceToAlign.txt

muscle -align sequenceToAlign.txt -output msa_align.fa &> outMuscle 
