#!/bin/bash

echo "Parse trimmed ends"
cd data
perl ~ssajjadi/wgacbin/step_8_mpi/collectTrim.pl step_8_mpi/trim/
mv step_8_mpi/trim/ParallelOutput.trim both.parse.defugu.trim
rm -rf step_8_mpi

echo "Fix overlapping alignments"
blast_align_by_align_overlap_fix3.pl -i both.parse.defugu.trim

echo "Defractionate alignments"
blast_defractionate3.pl -s 400000 -t both.parse.defugu.trim.fixed.trim -f ../fastawhole

echo "Now run global alignments on eeek:"
echo "  /net/eichler/vol5/home/tinlouie/wgacscripts/wgacAlign.sh ../fasta 64bit"