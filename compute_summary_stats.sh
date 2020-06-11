module load bedtools/2.29.0

cat GenomicSuperDup.tab | awk '{if ($3-$2 >= 1000 ) print }' | awk '{if ( $26 >= 0.9 ) print $0 }' >temp1000_90percent.tab
grep "^chr" temp1000_90percent.tab | awk '{ if ( $4 ~ /^chr/ ) print }'  >temp1000_90percent_just_chr.tab

echo "number of bases in all of these seg dups (not counting any bases twice)"
bedBases.sh <temp1000_90percent.tab
echo "just chr:"
bedBases.sh <temp1000_90percent_just_chr.tab

echo "number of non-redundant loci:"
cat temp1000_90percent.tab | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
echo "just chr (regardless if the other end is in unplaced):"
cat temp1000_90percent_just_chr.tab | awk '{print $1"\t"$2"\t"$3}' | grep -v "ML" | grep -v "QNV" | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l


echo "number of pairwise alignment 1kb and 90% and 1kb:"
GENOMIC=`wc -l temp1000_90percent.tab | awk '{print \$1}' `
echo $(( $GENOMIC / 2 ))
echo "just chr:"
GENOMIC_JUST_CHR=`wc -l temp1000_90percent_just_chr.tab | awk '{print \$1}' `
echo $(( $GENOMIC_JUST_CHR / 2 ))


# interspersed alignments

echo "number of pairwise alignments between different chromosomes"
DIFF_CHROMOSOMES=`cat temp1000_90percent_just_chr.tab | awk '{if ( \$1 != \$7 ) print }' | wc -l | awk '{print \$1}'`
echo $(( $DIFF_CHROMOSOMES / 2 ))

echo "number of pairwise alignments on same chromosome but at least 1 MB apart:"
SAME_CHROMOSOME_1MB_APART=`cat temp1000_90percent_just_chr.tab | awk '{if ( \$1 == \$7 ) print }' | awk '{ if ( \$2 < \$8 ) nDist = \$8 - \$3; else nDist = \$2 - \$9; if ( nDist >= 1000000 ) print }' | wc -l | awk '{print $1}' `
echo $(( $SAME_CHROMOSOME_1MB_APART / 2 ))

echo "interspersed seg dups:"
echo $(( $DIFF_CHROMOSOMES / 2 + $SAME_CHROMOSOME_1MB_APART / 2 ))



