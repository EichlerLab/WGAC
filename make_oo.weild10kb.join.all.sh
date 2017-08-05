#$ -S /bin/sh
#!/bin/bash

module load modules modules-init modules-gs

#usage qsub ./blastparser42.sh `pwd`/blastout `pwd`/data/blastout.parse 250

INPUT="$1"
OUTPUT="$2"

in_name=`echo "$INPUT"|cut -d "/" -f2`

cd data

perl /net/eichler/vol2/local/inhousebin/blast_hit_by_hit_welder2.pl -i $in_name -g 40:10000:20 -o oo.weild10kb 

cd ../
