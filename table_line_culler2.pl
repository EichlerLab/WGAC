#!/usr/bin/perl
#LOAD MODULES

use Getopt::Std;
use strict 'vars';

use vars qw($program $pversion $pdescription $pgenerate);
$program = "$0";
$program =~ s/^.*\///;

### program stats ###
$pversion ='2.001019';
$pdescription="$program (ver:$pversion) "
. "allows for the culling of certain rows from a table based on "
. "any criteria designated.  It can also be used to modify rows or columns. "
. "In fact it is pratically a spreadsheet except for the graphical interface and "
. "ease of use.";
$pgenerate= 'jeff:dnhc genetics:dnh';
### program stats end ###
#print "usage: $program -in [path] -out [path] [options]\n";


####000814 added program,version, description 
#### changed output 0 or 1#####
#### added hashtable variable######


use vars qw($true $false);
($true,$false)=(1,0);
use vars qw($opt_p $opt_h $opt_k $opt_o $opt_t $opt_v $path @files);
use vars qw(%h @c);

if (!defined $ARGV[0]) {

print "USAGE\n$program -p [file or dir path] -t '[test critera]' [options]
DESCRIPTION\n$pdescription
REQUIRED INPUT
-p [path] to file or directory of files to process
-t the test critera in quotes(' ') to be evaluated for truth
   *use \$c[1] for column variable (columns begin at 0)
   *the critera is tested for trueness 0 is false any other number true
   *\$c[1]>100  <====column 1 must be greater than 100 for inclusion
   *\$c[1]+=100; return 1; <====column 1 + 100 nothing culled since returns true;
OPTIONAL INPUT
-h [path] to file of tab-delimited data to place in hash
   *use \$h{key}[0] to acess column zero of a particular hash key line
  -k [number] of column in table from which to build hash (default is 0)	 
-o [directory] in which to place out files (default same as files)
-v [switch] verbose for debugging or to appear hard-working
   *prints 'T' and 'blanks' so that you can see the pattern of culling
";
exit;
}
if ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
	system "perldoc $0\n";
	exit;
}


getopts('p:t:o:h:k:v');
defined $opt_p || die "Please use -p to input path to file or directory\n";
defined $opt_t || die "Please use option -t to designate test critera\n";

#$opt_t=' $c[1]+=$h{$c[0]}[1];$c[2]+=$h{$c[0]}[1];$c[5]+=$h{$c[4]}[1];$c[6]+=$h{$c[4]}[1]; !($h{$c[0]}[1] && $h{$c[4]}[1])' ;

if (opendir (DIR, $opt_p) ) {
	@files = grep { /[a-zA-Z0-9]/ } readdir DIR;
	close DIR;
	$path = $opt_p;
} elsif (open (IN, $opt_p)) {
	($path ) = $opt_p =~ /(^.*)\//;
	$opt_p=~ s/^.\///;
	$path ||='.';
	@files=($opt_p);
} else {
	die "-p  ($opt_p)  not a file and not a directory\n";
}
defined $opt_o || ($opt_o = $path);
####load a hashtable####
%h=();
if ($opt_h) {
	print "HASHTABLE LOADING\n";
	open(HASHTABLE, $opt_h) || die "Can't open hast table ($opt_h) for reading!\n";
	my $header=<HASHTABLE>;
	while (<HASHTABLE>) {
		s/\r\n/\n/;
		chomp;
		next if !/\w/;
		my @c=split "\t";
		print "$c[$opt_k]\n";
		$h{$c[$opt_k]}=\@c;
	}
	print "HASHTABLE LOADED (",scalar(keys %h)," KEYS)\n";
}
#####cull -p files#####
foreach my $f (@files) {
	open (OUT, ">$opt_o/$f.cull") || die "Can't opent outfile ($opt_o/$f.tmp)!\n";
	print "******$opt_o/$f" if $opt_v;
	open (FILE,"$path/$f") || die "Can not open file $path/$f\n";
	my $header=<FILE>;
	print OUT $header;
	my $evaluated;
	my @c;
	my $line;
	while (<FILE>) {	
      s/\r\n/\n/; chomp; chomp;
		@c=split "\t";
		$line=$_;
		$evaluated= eval ($opt_t);
		if ($evaluated) {
			print "T" if $opt_v;
			print OUT (join("\t",@c),"\n");
		}else {
			print " " if $opt_v;
		};
		
	}
	close FILE;
	close OUT;
}




