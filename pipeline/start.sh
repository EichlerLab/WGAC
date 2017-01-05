#!/bin/bash

if [[ "$#" -lt "2" ]]
then
  echo "Usage: start.sh <species> <chrom_dir>"
  exit 1
fi

SPECIES=$1
CHROMOSOMES_DIR=$2

mkdir -p wgac
cd wgac

ln -s $CHROMOSOMES_DIR fastawhole

echo "Fractionate"
mkdir -p fasta
fasta_fractionate -f fastawhole -s 400000 -o fasta

echo "Get mask out from lowercase letters"
mkdir -p mask_out
ls fasta | xargs -i perl /net/eichler/vol4/home/jlhudd/wgac/maskOutGenFromLowCase.pl fasta/{} mask_out/{}.out

echo "Fuguize sequences"
mkdir -p fugu
fasta_fuguize_batch.pl -f fasta -r mask_out -o fugu

echo "Create blastdb"
mkdir -p blastdb
dir_cat fugu blastdb/bofugu
formatdb -i blastdb/bofugu -o F -p F -a F
rm -f blastdb/bofugu 

echo "Submit blast job to cluster"
mkdir -p blastout
qsub -N "${SPECIES}_wgac" -q all.q /net/eichler/vol4/home/jlhudd/wgac/blast64.sh -t "/var/tmp/jlhudd/${SPECIES}_wgac"

echo "Prepare fugu sequences for self-blast"
cp -r fugu fugu2
fasta_only_ATCGN.pl -p fugu2 -o fugu2 -u
mkdir -p selfblast

echo "Now start self-blast on srna3.gs.washington.edu:"
echo "  ssh srna3.gs.washington.edu"
echo "  perl fugumation.pl -i `pwd`/fugu2 -o `pwd`/selfblast"