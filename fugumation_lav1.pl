#!/usr/bin/perl
#####THIS LITTLE SCRIPT MUST BE RUN WITH IN THE DIRECTORY CONTAINING THE DIRS:
#########  FASTA, MASK_OUT;
use strict 'vars';
use Getopt::Std;


use vars qw($opt_m $opt_h $opt_l);
use vars qw($true $false);
($true,$false)=(1,0);

getopts('mhl');
if ( $opt_h) {
print "usage: .pl [options] **************
OPTIONAL PARAMETERS
-h  help (display this text)
-m  use multiple processors (two to be precise)
-l  skip execution of webb_batch
	
";
exit;
}



print "OPERATING SYSTEM: $^O\n";  ##THIS MUST RUN ON UNIX##

#################################
########   batch_start.pl    must be run first!!!!
#################################
#################################

system "mkdir data";



######################################################
#######LAV ANALYSIS  ###################################
###################################################
die "No fugu directory! Have you run batch_start.pl yet?" if ! -e "fugu";

#get simple character set#
if (!$opt_l) {
	system "fasta_only_ATCGN.pl -p fugu -o fugu -u";


	####THIS PROVIDES THE COMMAND IN A FILE SO MISTAKES WON'T BE MADE###
	open (TEMP, ">lavcom.tmp") || die "Can't open selfcom.tmp\n!";
	print TEMP "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-75 V=-75 O=180 E=1 W=14 Y=1400'";
	close TEMP;


	system "rm -rf lav_int";
	system "mkdir lav_int";
	#attempt 1#
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-75 V=-75 O=180 E=1 W=14 Y=1400' ";
	#attempt 2#
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-75 V=-75 O=180 E=1 W=14 Y=1400'";
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-75 V=-75 O=180 E=1 W=14 Y=1200'";
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-75 V=-75 O=180 E=1 W=14 Y=1000'";
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-75 V=-75 O=180 E=1 W=14 Y=900'";
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-75 V=-75 O=180 E=1 W=14 Y=800'";
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-75 V=-75 O=180 E=1 W=14 Y=600'";
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-80 V=-80 O=180 E=1 W=14 Y=400'";
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-80 V=-80 O=180 E=1 W=14 Y=200'";
	system "webb_batch.pl -f fugu -o lav_int -p 'B=2 M=30 I=-90 V=-90 O=180 E=1 W=18 Y=200'";
}
#I didn't need last one!

#my $pause=<STDIN>;



######BREAK OVERLAPPING LAVS #######
print "BREAKING SELF OVERLAPS UP\n";
system "mkdir lav_int2";
system "blast_lav_break_self_overlap2.pl --in lav_int --out lav_int2";



#####PARSE THE WONDERFUL WORLD OF LAV ########
print "PARSING...\n";
system "blast_lav_hit_by_hit.pl --in lav_int2 --out data/lav_int2.parse "
		." -options 'MIN_BPALIGN=>200, MIN_FRACBPMATCH=>0.88, MAX_%GAP => 40, SKIP_SELF => 0, SKIP_OVERLAP=>1'";









