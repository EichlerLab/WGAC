#!/bin/bash
echo "Are you on a node with >2G of RAM? (y/n)"
read confirmation

if [[ "$confirmation" = "n" ]]
then
    echo "Ok, get a node with some more memory."
    exit 1
fi

echo "Parse trimmed ends"
cd data

if [[ -e step_8_mpi/trim/ ]]
then
    perl ~ssajjadi/wgacbin/step_8_mpi/collectTrim.pl step_8_mpi/trim/
fi

if [[ -e step_8_mpi/trim/ParallelOutput.trim ]]
then
    mv step_8_mpi/trim/ParallelOutput.trim both.parse.defugu.trim
    rm -rf step_8_mpi
fi

echo "Fix overlapping alignments"
blast_align_by_align_overlap_fix3.pl -i both.parse.defugu.trim

echo "Defractionate alignments"
blast_defractionate3.pl -s 400000 -t both.parse.defugu.trim.fixed.trim -f ../fastawhole

echo "Now run global alignments on eeek:"
echo "  /net/eichler/vol5/home/tinlouie/wgacscripts/wgacAlign.sh ../fasta 64bit"