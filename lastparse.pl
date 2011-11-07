#!/usr/bin/perl
#####THIS LITTLE SCRIPT MUST BE RUN WITH IN THE DIRECTORY CONTAINING THE DIRS:
#########  FASTA, MASK_OUT;


#####PARSE THE WONDERFUL WORLD OF LAV ########
print "PARSING...\n";
system "/home/linchen2/bin/wgac/blast_lav_hit_by_hit.pl --in lav_int2 --out data/lav_int2.parse "
		." -options 'MIN_BPALIGN=>200, MIN_FRACBPMATCH=>0.88, MAX_%GAP => 40, SKIP_SELF => 0, SKIP_OVERLAP=>1'";









