#!/bin/bash


if [ $# -ne 3 ]
then
 echo "please give args 1.filename, 2.dbname, 3.outfile  name"
 exit
fi

loadOOWELD2Mysqltable.pl $1 $2 tmp

mysql -plinchen2 -ulinchen2 $2 -e"
select QNAME as chrom, QB as chromStart, QE as chromEnd, concat(sname, ':', least(sb, se)) as name, score, if(sb<se, '+', '_') as strand, sname as otherChrom, least(sb, se) as otherStart, greatest(sb,se) as ohterEnd, slen as otherSize, 0 as uid, 1000 as posBasesHit, 'N/A' as testResult, 'N/A' as verdict, 'N/A' as chits, 'N/A' as ccov, file as alignfile, base_s as alignL, indel_n, indel_s, bpalign as alignB, base_match as matchB, base_mis as mismatchB, transitions, transversions, per_sim as fracMatch, per_sim_indel as fracMatchIndel, k_jc as jcK, k_kimura as k2K from tmp
union
select SNAME as chrom, least(SB,SE) as chromStart, greatest(SB,SE) as chromEnd, concat(qname, ':', qb) as name, score, if(sb<se, '+', '_') as strand, qname as otherChrom, qb as otherStart, qe as ohterEnd, qlen as otherSize, 0 as uid, 1000 as posBasesHit, 'N/A' as testResult, 'N/A' as verdict, 'N/A' as chits, 'N/A' as ccov, file as alignfile, base_s as alignL, indel_n, indel_s, bpalign as alignB, base_match as matchB, base_mis as mismatchB, transitions, transversions, per_sim as fracMatch, per_sim_indel as fracMatchIndel, k_jc as jcK, k_kimura as k2K
from tmp" > $3genomicSuperDup.tab
rm_header.pl $3genomicSuperDup.tab

