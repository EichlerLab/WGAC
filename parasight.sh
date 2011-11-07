#!/bin/bash
#
# Uses miropeats blast output and genomic super dups to find overlaps between
# BACs and WGAC positive regions. Output is a Parasight extra file which stacks
# WGAC positive regions by chromosome.

# Create a BED file with chromosome name, start, and end as the first three
# columns and the BAC values as the next three. The final column is the %
# identity. Alignment length has already been filtered by the miropeats pipeline
# at this point.
awk 'OFS="\t" { if ($9 >= 0.94) { print $5,$6,$7,$1,$2,$3,$9 } }' parse_901k_sort.mod | sed 1d > blasthits.bed

# Run intersectBed between blasthits and the WGAC superdup output which includes
# chrom and otherChrom values. Only keep the otherChrom and the BAC values.
echo Enter location of super dup file:
read superdup
intersectBed -a blasthits.bed -b $superdup -wb | cut -f 4,5,6,14 | sort -k 1,1 -k 2,2n -k 3,3n -k 4.4,4.5n > bac.matches

# Split output into one file per chromosome, merge the coordinates per file, and
# rejoin the merged files.
#
# Also sort the merged BAC output by BAC name, chromosome, BAC start, and BAC
# end.
for chrom in `cut -f 4 bac.matches | sort | uniq`
do
    grep -w $chrom bac.matches | mergeBed -i stdin -d 1 | awk -v chrom="$chrom" 'OFS="\t"{ print $0,chrom }'
done | sort -k1,1 -k4,4 -k2,2n -k 3,3n > bac.matches.merged

# Filter BAC matches extra file to keep only contigs whose BAC/chromosome pair
# has at least one segment >= 5kb.
awk 'OFS="\t"{if ($3 - $2 >= 5000 || $3 - $2 <= -5000) {print $1, $4}}' bac.matches.merged | sort | uniq > seq_chr.uniq

# Create a list of unique chromosomes and assign a unique color to each one in a
# file called "chromosomes_color".
cut -f 2 seq_chr.uniq | sort | uniq > chromosomes
python ~/parasight/colorgen.py

# Create a Parasight extra file for the matches.
python ~/parasight/writeExtraFile.py --input=bac.matches.merged \
--output=bac.matches.extra --color=map --width=4 --type=WGAC \
--extraheaders="feature" --feature="feature"

# Stack the matching WGAC results by chromosome to prevent overlapping segments.
python ~/parasight/stack_extras.py --filter seq_chr.uniq \
  bac.matches.extra bac.stacked.extra \
  --featurename="feature" \
  --offset=22 --increment=10