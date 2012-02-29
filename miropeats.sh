#!/bin/bash

if [[ "$#" -lt "3" ]]
then
    echo "Usage: ~/wgac/miropeats.sh full_species wgac_superdup path_to_chromosomes"
    exit 1
fi

FULL_SPECIES=$1
SUPERDUP=$2
CHROMOSOMES_DIR=$3

# Get all overlapping regions between BACs and WGAC keeping only alignment
# data (like a BLAST report).
#
# Cut 4-6,8-10 instead here to get "true mapped" regions.
intersectBed -a bac_coordinates.bed -b $SUPERDUP -wb | \
  cut -f 4-6,14-16 | sort | uniq > bac_wgac.unfiltered.alignment
 
# Get valid (size-filtered) BAC/chromosome regions for WGAC.
sed 1d bac.matches.extra | cut -f 1-4 | sort | uniq > bac_wgac_valid_pairs
 
# Get chromosome coordinates for valid BAC/chromosome pairs by intersecting
# on BAC coordinates and only keep matches where chromosomes (columns 4 and 8)
# are the same.
intersectBed -a bac_wgac_valid_pairs -b bac_wgac.unfiltered.alignment -wa -wb | \
  awk '$4 == $8' | cut -f 1-3,8-10 > bac_wgac.filtered.alignment
 
# Get just BAC names and chromosome regions.
cut -f 1,4-6 bac_wgac.filtered.alignment | sort | uniq > \
  bac_wgac.filtered.regions
 
# Merge regions across 250 Kbp for each BAC.
for bac in `cut -f 1 bac_wgac.filtered.regions | sort | uniq`
do
  grep "$bac" bac_wgac.filtered.regions | cut -f 2-4 | sort | uniq | \
    mergeBed -i stdin -d 250000 | \
    awk -v bac="$bac" 'OFS="\t"{ print bac,$0 }'
done > bac_wgac.merged.regions

# Get chromosome regions coordinates.
cut -f 2-4 bac_wgac.merged.regions | sort | uniq > chromosomes.bed
 
# Add ".fa" to the first column's values to make fasta_subseq work
# with .fa files.
sed 's/\(chr\S\+\)/\1.fa/' chromosomes.bed > chromosomes.fa.bed

# Cut regions.
mkdir -p chromosomes
fasta_subseq -L chromosomes.fa.bed -h \
  -D $CHROMOSOMES_DIR \
  -o chromosomes/
 
# Remove the .fa from each filename.
rename '.fa_' '_' chromosomes/chr*

# Repeatmask chromosomes.
cd chromosomes
ls chr* > repeats.txt
RepeatMasker -no_is -xsmall -qq -species "$FULL_SPECIES" -dir `pwd` `cat repeats.txt`
 
# Clean up.
rm -f *.ref *.tbl *.cat *.log *.masked
cd ..
