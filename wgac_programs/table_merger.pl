#!/usr/bin/perl

use Getopt::Std;
use strict 'vars';

use vars qw($true $false);  ($true,$false)=(1,0);

use vars qw($opt_i $opt_s $opt_o);
#use vars qw();

if (!defined $ARGV[0]) {
print "USAGE table_merger.pl -i [TABLE1]:[TABLE2] -s '[search (parm1)]:[search (parm2)]'
This program merges TABLE2 to TABLE1.  If TABLE1 is missing TABLE2 row TABLE2 is lost.
REQUIRED OPTIONS:
-i [path] to file or directory of files to process
-s the test critera  in quotes (' ') to be evalutated use () for actual part to match
-o output table path
";
exit;
}

getopts('i:s:o:');

$opt_i =~/\:/ || die "Please use -i with [PATH TABLE1]:[PATH TABLE 2].\n";
$opt_s =~/\:/ || die "Please use -s with 'pattern1:pattern2' ($opt_s)\n";
defined $opt_o || die "Please use -o to enter path of resulting merged table.\n";


my ($table1,$table2) = split ':', $opt_i;
my ($pattern1,$pattern2) = split ":",$opt_s;
print "$pattern1---$pattern2 \n";
#my $pause=<STDIN>;
open (TABLE2, "$table2") || die "Can't open table2 for reading ($table2)!\n";
my $header2 = <TABLE2>;
my %h2;
my $line=0;
print "LOADING $table2 AS HASH!\n";
while (<TABLE2>) {
	chomp;
	#print $_, "\n";
	if ( /$pattern2/ ) {
		my $x=$1;
		$x =~s/^0+//;
		print "HERE IT IS:$x \n";
		$h2{$x}=$_;
		
	}
}
#my $pause=<STDIN>;
open (OUT, ">$opt_o") || die "Can't open output table ($opt_o)!\n";

open (TABLE1, "$table1") || die "Can't open table1 for reading ($table1)!\n";
my $header1 = <TABLE1>;
chomp $header1;
print OUT "$header1\t$header2";

print "JOINING HASH TO $table1\n";
while (<TABLE1>) {
	chomp;
	print OUT;
	if ( /$pattern1/ ) {
		my $x=$1;
		$x=~s/^0+//;
		print"TABLE1:$x\n";
		
		print OUT "\t$h2{$x}" if defined $h2{$x}
	}
	print OUT "\n";
}
close OUT;

