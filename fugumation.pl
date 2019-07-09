use strict 'vars';
use Getopt::Std;

# usage:
# perl fugumation.pl -i fugu -o selfblast
# 	fugu is name of directory containing fugu files
# 	selfblast is name of directory to hold output files

use vars qw($opt_i $opt_o);
getopts('i:o:');

#die "cannot create subdir lav_int \n" unless (mkdir "$opt_o/lav_int");
#die "cannot create subdir lav_int2 \n" unless (mkdir "$opt_o/lav_int2");
#die "cannot create subdir data \n" unless (mkdir "$opt_o/data");


# replaces RYKMSWBDHV with N in a fasta record, switch to uppercase
#system "perl fasta_only_ATCGN.pl -p $opt_i -o $opt_i -u";

runWebb(-75, 14, 1400);
runWebb(-75, 14, 1200);
runWebb(-75, 14, 1000);
runWebb(-75, 14, 900);
runWebb(-75, 14, 800);
runWebb(-75, 14, 600);
runWebb(-80, 14, 400);
runWebb(-80, 14, 200);
runWebb(-90, 18, 200);

######BREAK OVERLAPPING LAVS #######
# system "perl blast_lav_break_self_overlap2.pl --in $opt_o/lav_int --out $opt_o/lav_int2";

#####PARSE THE WONDERFUL WORLD OF LAV ########
# system "perl blast_lav_hit_by_hit.pl --in $opt_o/lav_int2 --out $opt_o/data/lav_int2.parse -options 'MIN_BPALIGN=>200, MIN_FRACBPMATCH=>0.88, MAX_%GAP => 40, SKIP_SELF => 0, SKIP_OVERLAP=>1'";

exit;

# webb_batch.pl was boiled down to the following subroutine:
# invokes the webb_self program on each file in fugu directory, hopefully creates corresponding output file
# ... won't create output file if it already exists (from a previous iteration)
# ... removes 'core' file which is the result of segmentation fault
sub runWebb {
	my $IVparameter = shift; # penalty
	my $Wparameter = shift; # word size
	my $Yparameter = shift; 

	opendir (DIR, $opt_i) || die "cannot open directory to fugu files\n";

	while ( defined (my $f = readdir(DIR)) ) {
	    next if $f =~ /^\./;
		my $outputFile = "$opt_o/$f.intf";
		
		next if (-e $outputFile); # OUTPUT FILE ALREADY EXISTS
		
		print "DOING $f...\n";
		system "./webb_self $opt_i/$f B=2 M=30 I=$IVparameter V=$IVparameter O=180 E=1 W=$Wparameter Y=$Yparameter > $outputFile";

		# if output file exists but its first line is "bad", then delete the file
		# so that it can be redone in the next round
		if (open (POSTCHECK, $outputFile)) {
			my $line= <POSTCHECK>;
			close POSTCHECK;
			unlink $outputFile unless ( $line =~ /\#/ );
		}

		# if webb_self had a core dump, remove the file (it can grow big)
		unlink 'core' if (-e 'core');
		
	}

	closedir DIR;
	
}
