cat fastalength.log | sed 1d | sed '$d' | grep -vw -f chromosomes.txt | sort -k3,3nr | awk 'BEGIN {x = 0} {print $2, $3, x; x += $3}' >fakeChromosomeOffsets.txt

