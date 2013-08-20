#!/bin/bash

if [ -z $1 ]
then
  echo "Usage: nr_report.sh fastalength.log interDup intraDup AllDupLen <output.tab> [factor]"
  echo
  echo "    factor: the number of base pairs to divide each duplication value by. Default: 1000"
  echo "    For example: factor=1000 for duplication in Kbp, factor=1000000 for Mbp"
  exit 1;
fi

FASTALENGTH=$1
INTER=$2
INTRA=$3
ALL=$4
OUTPUT=$5

if [ -z $6 ]
then
    FACTOR=1000
else
    FACTOR=$6
fi

# Get a list of all contigs.
cut -f 2 "$FASTALENGTH" | sed '1d;$d' | sort -k 1,1 > "$OUTPUT"
sort -k 1,1 -o ${INTER} ${INTER}
sort -k 1,1 -o ${INTRA} ${INTRA}
sort -k 1,1 -o ${ALL} ${ALL}

# Join all files on contig name.
join -t '	' -a 1 -j 1 -o 1.1 2.2 -e "0" "$OUTPUT" "$INTER" > tmp && mv -f tmp "$OUTPUT"
join -t '	' -a 1 -j 1 -o 1.1 1.2 2.2 -e "0" "$OUTPUT" "$INTRA" > tmp && mv -f tmp "$OUTPUT"
join -t '	' -a 1 -j 1 -o 1.1 1.2 1.3 2.2 -e "0" "$OUTPUT" "$ALL" > tmp && mv -f tmp "$OUTPUT"

awk -v factor="$FACTOR" \
'BEGIN { OFS="\t"; print "contig","inter","intra","total" }
       { print $1,$2/factor,$3/factor,$4/factor }
' "$OUTPUT" > tmp && mv -f tmp "$OUTPUT"