#$ -S /bin/sh
#!/bin/bash

module load modules modules-init modules-gs

#usage qsub ./blastparser42.sh `pwd`/blastout `pwd`/data/blastout.parse 250

INPUT="$1"
OUTPUT="$2"

cd data
table_line_culler2.pl -p oo.weild10kb.pieces.all -t'$c[27]>=0.90 && $c[64]>=1000'
table_line_culler2.pl -p oo.weild10kb.join.all -t'$c[27]>=0.90 && $c[22]>=1000'
cd ../
