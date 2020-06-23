#!/bin/bash

#script to estimate sequencing error rate
#input is a bamfile of reads mapped to a reference
#and the reference, a mitochondrial genome, assumes no in-species variation
#usage:
  #errorRate.sh <aln.bam> <ref.fa> <indels considered? 1/0 yes or no>

aln=$1
ref=$2

len=$(samtools view -h $aln | tail -n 1 |  cut -f10 | wc -c)

if [ $3 == 0 ]
then
  samtools view -h $aln | \
  awk '$1 ~ "^@" || $6 !~ "I|D"' | \
  samtools view -b - | \
  samtools calmd -e - $ref | \
  cut -f10 | \
  sed 's/=//g' > temp
else
  samtools calmd -e $aln $ref | \
  cut -f10 | \
  sed 's/=//g' > temp
fi

while IFS= read -r line; do echo ${#line}; done < temp > count

awk '{print $1 / '$len'}' count > $aln.error

rm count
rm temp
