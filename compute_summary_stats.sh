module load bedtools/2.29.0

cat GenomicSuperDup.tab | awk '{if ($3-$2 >= 1000 ) print }' | awk '{if ( $26 >= 0.9 ) print $0 }' >temp1kb_90percent.tab

cat temp1kb_90percent.tab | ../filterByTokenValue.py --szFileOfLegalValues ../chromosomes.txt --n0BasedToken 0 | ../filterByTokenValue.py --szFileOfLegalValues ../chromosomes.txt --n0BasedToken 6 >temp1kb_90percent_just_chr.tab

# above replaces this:
#grep "^chr" temp1kb_90percent.tab | awk '{ if ( $7 ~ /^chr/ ) print }' | grep -v random | grep -v chrUn  >temp1kb_90percent_just_chr.tab




cat GenomicSuperDup.tab | awk '{if ($3-$2 >= 1000 ) print }' | awk '{if ( $26 >= 0.95 ) print $0 }' >temp1kb_95percent.tab

cat temp1kb_95percent.tab | ../filterByTokenValue.py --szFileOfLegalValues ../chromosomes.txt --n0BasedToken 0 | ../filterByTokenValue.py --szFileOfLegalValues ../chromosomes.txt --n0BasedToken 6 >temp1kb_95percent_just_chr.tab


# above replaces this:
#grep "^chr" temp1kb_95percent.tab | awk '{ if ( $7 ~ /^chr/ ) print }' | grep -v random | grep -v chrUn >temp1kb_95percent_just_chr.tab


cat GenomicSuperDup.tab | awk '{if ($3-$2 >= 1000 ) print }' | awk '{if ( $26 >= 0.98 ) print $0 }' >temp1kb_98percent.tab


cat temp1kb_98percent.tab | ../filterByTokenValue.py --szFileOfLegalValues ../chromosomes.txt --n0BasedToken 0 | ../filterByTokenValue.py --szFileOfLegalValues ../chromosomes.txt --n0BasedToken 6 >temp1kb_98percent_just_chr.tab

# above replaces this:
#grep "^chr" temp1kb_98percent.tab | awk '{ if ( $7 ~ /^chr/ ) print }' | grep -v random | grep -v chrUn >temp1kb_98percent_just_chr.tab

for szFile in temp1kb_90percent.tab temp1kb_95percent.tab temp1kb_98percent.tab temp1kb_90percent_just_chr.tab temp1kb_95percent_just_chr.tab temp1kb_98percent_just_chr.tab
do
    echo $szFile

    echo "# of pairwise aligments:"
    GENOMIC=`wc -l $szFile | awk '{print \$1}' `
    echo $(( $GENOMIC / 2 ))

    echo "  inter:"
    GENOMIC=`cat $szFile | awk '{if ( \$1 != \$7 ) print }' | wc -l | awk '{print \$1}' `
    echo $(( $GENOMIC / 2 ))
    DIFF_CHROMOSOMES=$(( $GENOMIC / 2 ))

    echo "  intra:"
    GENOMIC=`cat $szFile | awk '{if ( \$1 == \$7 ) print }' | wc -l | awk '{print \$1}' `
    echo $(( $GENOMIC / 2 ))

    echo "number of pairwise alignments on same chromosome but at least 1 MB apart:"
    SAME_CHROMOSOME_1MB_APART2=`cat $szFile | awk '{if ( \$1 == \$7 ) print }' | awk '{ if ( \$2 < \$8 ) nDist = \$8 - \$3; else nDist = \$2 - \$9; if ( nDist >= 1000000 ) print }' | wc -l | awk '{print $1}' `
    echo $(( $SAME_CHROMOSOME_1MB_APART2 / 2 ))
    SAME_CHROMOSOME_1MB_APART=$(( $SAME_CHROMOSOME_1MB_APART2 / 2 ))

    # interspersed alignments

    echo "interspersed seg dups:"
    echo $(( $DIFF_CHROMOSOMES + $SAME_CHROMOSOME_1MB_APART ))

done


echo ""
echo "now bp"
echo ""



for szFile in temp1kb_90percent.tab temp1kb_95percent.tab temp1kb_98percent.tab temp1kb_90percent_just_chr.tab temp1kb_95percent_just_chr.tab temp1kb_98percent_just_chr.tab
do
    echo $szFile

    echo "# of bases:"
    cat $szFile | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{x += ( $3 - $2); y+=1 } END { print x }' | xargs printf "%'d\n"

    echo "  inter:"
    cat $szFile | awk '{if ( $1 != $7 ) print }' | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{x += ( $3 - $2); y+=1 } END { print x }' | xargs printf "%'d\n"


    echo "  intra:"
    cat $szFile | awk '{if ( $1 == $7 ) print }' | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{x += ( $3 - $2); y+=1 } END { print x }' | xargs printf "%'d\n"

done


echo ""
echo "now nonredundant loci"
echo ""



for szFile in temp1kb_90percent.tab temp1kb_95percent.tab temp1kb_98percent.tab temp1kb_90percent_just_chr.tab temp1kb_95percent_just_chr.tab temp1kb_98percent_just_chr.tab
do
    echo $szFile

    LOCI=`cat $szFile | awk '{print \$1"\t"\$2"\t"\$3 }' | sort -k1,1V -k2,2n | bedtools merge | wc -l | awk '{print \$1}'`
    echo $LOCI

done

echo ""
echo "inter-chromosomal"
echo ""


for szFile in temp1kb_90percent_just_chr.tab temp1kb_95percent_just_chr.tab temp1kb_98percent_just_chr.tab
do
    echo $szFile

    LOCI=`cat $szFile | awk '{if ( \$1 != \$7 ) print }' | awk '{print \$1"\t"\$2"\t"\$3 }' | sort -k1,1V -k2,2n | bedtools merge | wc -l | awk '{print \$1}'`
    echo $LOCI

done



echo ""
echo "intra-chromosomal"
echo ""


for szFile in temp1kb_90percent_just_chr.tab temp1kb_95percent_just_chr.tab temp1kb_98percent_just_chr.tab
do
    echo $szFile

    LOCI=`cat $szFile | awk '{if ( \$1 == \$7 ) print }' | awk '{print \$1"\t"\$2"\t"\$3 }' | sort -k1,1V -k2,2n | bedtools merge | wc -l | awk '{print \$1}'`
    echo $LOCI

done

