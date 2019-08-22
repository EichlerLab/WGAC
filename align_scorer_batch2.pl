#!/usr/bin/perl 

use strict 'vars';
use Getopt::Std;

use vars qw($default_executable);
$default_executable='align_scorer2.pl';


use vars qw($opt_a $opt_p $opt_o $opt_d $opt_e);
use vars qw( $temp  @dir);

if ($ARGV[0] eq "") {
	print "**************************************\n";
	print "-a  [directory] where alignment files are stored\n";
	print "-p  search pattern for files";;
	print "-o  name of summary output file\n";
	print "-d  [directory] where data should be stored\n";
	print "-e  [program path] chose a different program version\n";
	print "    (default is $default_executable)\n";
	print "\nYou must start in directory where you want output deposited!\n";
	die   "  ";
}

getopts('a:p:o:d:e:');
$opt_e ||= $default_executable;

system ("mkdir -p $opt_d") if $opt_d;

opendir (DIR, $opt_a) || die ('Cant open directory $opt_a!');
@dir = grep { /[A-Za-z0-9]/ } readdir (DIR);
close DIR;
@dir = sort @dir;
system "mkdir -p $opt_d";
unshift @dir, "";
print "DIR" ,join (":",@dir), "\n";
foreach my $d (@dir) {
	print "CHECKING $d\n";
	my $pathd="$opt_a/$d";
	next if ! (-d $pathd);
	print "  OPENING $d\n";
	opendir (SUBDIR,$pathd) || die ("Cant open directory $pathd!\n");
	my @filelist = grep { /$opt_p/i} readdir SUBDIR;
	close SUBDIR;
	print "  MAKE DIR:$d\n";
	system "mkdir $opt_d/$d";
	@filelist = sort @filelist;
	
	foreach my $file (@filelist) {
		system "$opt_e $pathd/$file";
		system "mv $pathd/$file.indel $opt_d/$d";
		system "mv $pathd/$file.mismatch $opt_d/$d";
	}
	#print "/nPAUSE";
	#my $pause=<STDIN>;
	
}

system "mv align_scorer.summary.table $opt_o";



