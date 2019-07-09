#!/usr/bin/perl
#00-06-01 renamed from fasta_encode_repeatmasker_batch.pl to fasta_fuguize_batch.pl


#!/usr/local/bin/perl
#LOAD MODULES

use Getopt::Std;
use strict 'vars';

use vars qw($true $false);
$true=1; $false=0;

use vars qw($opt_f $opt_r $opt_o $opt_g);
use vars qw(@fasta $ppath);
print "OPERATING SYSTEM: $^O\n";
$ppath="";
$ppath="perl /binp/" if $^O eq "MSWin32";
if ($ARGV[0] eq '') {

	print "fasta_fuguize_batch.pl ******************************************\n";
	print "This excutes a batch job of fasta_fuguize.pl\n";
	print "\t-f [path] directory containing fasta files\n";
	print "\t-r [path] directory containing repeatmasker *.out files\n";
	print "\t-o [path] direccotry to place encode fasta files *.encode\n";
	print "\t(the -g switch is alway on so proper nomenclature needed)\n";
	die "***********************************************************\n";
}

getopts('f:r:o:g');
$opt_f || die "Please enter with -f the path of a directory containing fasta files\n";
$opt_r || die "Please enter with -g the path for repeatmasker output\n";
$opt_o || die "Please enter with -o the output path for *.encode files\n";
opendir (DIR , $opt_f) || die "Can't open directory $opt_f\n";
@fasta =  grep {/[A-Za-z0-9]+/} readdir (DIR);
print @fasta, "\n";
close DIR;
my $count=0;
foreach (@fasta) {
	$count++;
	print "$count:$_ \n";
	my $f=$_;
	s/\.fa$//;
	s/\.fasta$//;
	system "$ppath"."fasta_fuguize.pl -f $opt_f/$f -r $opt_r/$_.out -o $opt_o/$_.fugu";
}
	



