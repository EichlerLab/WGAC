#!/bin/bash

echo "Clean up empty blast output"
perl /net/eichler/vol5/home/tinlouie/wgacscripts/rmNohitfiles.pl blastout

echo "Find incomplete blast results"
if find blastout -type f -exec grep -L 'Matrix' {} \;
then
  echo "Re-run blast to fix incomplete results"
  exit 1
fi

echo "Parse blast"
mkdir -p data
blastparser42.pl -in ./blastout \
  -out data/bo_self.parse \
  -noprevq -v \
  -output '-ssort=>name, -hsort=>qb' \
  -filter '-min_bpalign=>250, -min_fracbpmatch=>0.88, -max_%gap =>40, -no_subject_self => yes, -no_doubleoverlap=>score'

echo "Parse self-blast"
mkdir -p tmp
blast_lav_break_self_overlap2.pl --in selfblast --out tmp
mkdir -p data
blast_lav_hit_by_hit.pl --in tmp --out data/lav_int2.parse -options 'MIN_BPALIGN=>200, MIN_FRACBPMATCH=>0.88, MAX_%GAP => 40, SKIP_SELF => 0, SKIP_OVERLAP=>1'
rm -rf tmp
rm -rf fugu2

echo "Merge parsed blast results"
cd data
mv bo_self.parse both.parse
sed 1d lav_int2.parse >> both.parse
cd ..

echo "Defuguize"
blast_defuguize_hit_by_hit.pl -t data/both.parse -d .

echo "Prepare for end trimming"
cd data
mkdir -p step_8_mpi/defugu

# Second argument is number of lines per output file. The third argument is the
# number of lines in the defugu file.
DEFUGU_LINES=`wc -l both.parse.defugu | sed 's/\s\+/\t/g' | cut -f 1`
perl ~jlhudd/wgac/split.pl both.parse.defugu 300 $DEFUGU_LINES

find newdir/ -type f -exec cp {} step_8_mpi/defugu \;
rm -rf newdir
mkdir -p step_8_mpi/trim
cd ..

echo "Submit end trimming job to cluster"
qsub -N trim_ends -q all.q ~jlhudd/wgac/trim-ends.sh