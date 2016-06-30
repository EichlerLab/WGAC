#!/usr/bin/perl
# 00-05-21 modified to allow for case options -l -u and unchanged

use Getopt::Std;
use strict 'vars';

use vars qw($true $false);
use vars qw($opt_p $opt_l $opt_o $opt_v $opt_u $opt_l $opt_x);
use vars qw($path @files);
use vars qw($header @hheader %hash @c $badchar);

if ($ARGV[0] eq '') {

print "file_only_ATCGN.pl ******************************************
This replaces RYKMSWBDHV with N in a fasta record.
OPTIONS
  -p [path] to file or directory of files to process
  -o [directory] default is replace
  -u [switch] uppercase (default no case change)
  -l [switch] lowercase (default no case change)
  -v [switch] verbose for debugging
  -x [switch] do not remove extra characters
";
print "***********************************************************\n";
	exit;
}
getopts('p:o:vlux');
defined $opt_p || die "Please use -p to input path to file or directory\n";
$opt_o || ($opt_o = '.');


if (opendir (DIR, $opt_p) ) {
	@files = grep { /[a-zA-Z0-9]/ } readdir DIR;
	close DIR;
	$path = $opt_p;
} elsif (open (IN, $opt_p)) {
	($path ) = $opt_p =~ /(^.*)\//;
	$path||='.';
	$opt_p=~ s/^.\///;
	@files=($opt_p);
} else {
	die "-p  ($opt_p)  not a file and not a directory\n";
}


foreach my $f (@files) {
	open (OUT, ">$opt_o/$f.tmp") || die "Can't opent outfile ($opt_o/$f.tmp)!\n";
	print "$opt_o/$f.tmp" if $opt_v;
	open (FILE,"$path/$f") || die "Can not open file $path/$f\n";
	my $changes=0;
	while (<FILE>) {	
		s/\r\n/\n/;
		if ( /^>/ ) {
			print OUT $_;
			next;
		} else {
			$changes += tr/bdhvBDHVrykmRYKMswSW/nnnnNNNNnnnnNNNNnnNN/;
			$badchar+= tr/actgnACTGN//cd if !$opt_x;
			$badchar--;

		}
		if ($opt_u) {
			print OUT uc($_),"\n";
		} elsif ($opt_l) {
			print OUT lc($_),"\n";
		} else {
			print OUT $_,"\n";
		}
	}
	print "$path/$f => $changes ambigious changes,   $badchar bad characters removed\n";
	close FILE;
	close OUT;
	system "mv $opt_o/$f.tmp $opt_o/$f";
}




