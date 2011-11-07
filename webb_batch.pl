#!/usr/bin/perl
####THIS PROGRAM BLASTS A DIRECTORY OF FASTAS AMONGST THEMSELVES ###


#!/usr/local/bin/perl
#LOAD MODULES

use Getopt::Std;
use strict 'vars';

use vars qw($true $false);
use vars qw($program $pversion $pdescription);

$program = 'webb_batch';
$pversion = '0.001111';
$pdescription = "$program (ver:$pversion)".
' runs webb_self on a directory of fasta files with rudementary multiprocessor support\n';

$true=1; $false=0;

use vars qw($opt_f $opt_o $opt_p);
use vars qw($title $position $subseq $basename $path @fastas);

if ($ARGV[0] eq '') {

	print "$program ***************\n";
	print "DESCRIPTION\n$pdescription\n";
	print "ARGUMENTS\n";
	print "-f [path] directory of fasta files (use encode files)\n";
	print "-o [path] directory for output (*.intf)\n";
	print "-p [program parameters] enclosed in quotes e.g. for fuguization 'B=2 M=30 I=-75 V=-75 O=150 E=1 W=12 Y=1000'\n";
	die "***********************************************************\n";
}

getopts('f:o:p:');
$opt_f || die "Please input with -f  a directory of fasta files to concatenate and blast\n";
$opt_o || die "Please input with -o a basename for output\n";
$basename=$opt_o;
$basename=~s/(^.*)\///;
$true=1; $false=0;$true=1; $false=0;$true=1; $false=0;$path=$1;

print "RUNNING WEBB SELF ON DIRECTORY ($opt_f)...\n";
opendir (DIR, $opt_f) || die "Can't open directory $opt_f\n";
@fastas= sort grep { /[a-zA-Z0-9]/ } readdir (DIR);
closedir DIR;
open (DIROUT, $opt_o) || die "Out directory doesn't exist $opt_o\n";
closedir DIROUT;
my ($count);
foreach my $f (@fastas) {
  $count++;
  #print "$count  ";
  if (-e "$opt_o/$f.intf") {
	 print "$count) OUTPUT FILE ALREADY EXISTS..$f.skipping\n";
	 next;
  }
  print "$count) DOING $f...\n";
  system "webb_self $opt_f/$f $opt_p > $opt_o/$f.intf";
  #system "fasta_length2.pl $opt_f/$f";
  if  (open (POSTCHECK,"$opt_o/$f.intf")) {
	 my $line= <POSTCHECK>;
	 #print "    $line";
	 next if  $line =~ /\#/; #good output or at least the beginnings 
	 unlink "$opt_o/$f.intf"; #delete bad output and so it can be redone
  } 
  for (my $delay=1;$delay<100000; $delay++) { } #pause to allow a control c exit
}
for (my $delay=1;$delay<1000000; $delay++) { } #pause to allow another process to catch up
