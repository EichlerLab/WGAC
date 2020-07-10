sort -k1,1 -k2,2n | bedtools merge -i stdin | awk '{x += ( $3-$2); y+=1 } END { print x }' | xargs printf "%'d\n" 

