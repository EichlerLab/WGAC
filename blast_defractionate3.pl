#!/usr/bin/perl
####THIS PROGRAM UPDATES CONTIGS


#!/usr/local/bin/perl
#LOAD MODULES
#00-09-20  add option to simply work on 
#00-08-06  defractionate must beable to handle multiple directories####
#01-01
#00-05-26  fixed error in name parsing so that NT_num_(\d+) would capture just end digits
use Getopt::Std;
use strict 'vars';

use vars qw($true $false);
$true=1; $false=0;

use vars qw($opt_f $opt_s $opt_o $opt_t $opt_c $opt_l $opt_1 $opt_2);
use vars qw($title $position $subseq $header $count %len);

if ($ARGV[0] eq '') {

	print "blast_defractionate ******************************************\n";
	print "  This program decodes a pairwise blast table or single.\n";
	print "\t-t [path] a blast parse table file is required\n";
	print "\t-c [switch] no calculation of overall fasta length (-f not needed)\n";
	print "\t-f [dir] with whole fasta files\n";
	print "   (if diff dir for sub and query use [querypath:subpath]\n   Note:this is a quick fix shoud do better search\n";
	print "\t-s [integer] size of fractions to be reunited\n"; 
	print "\t-o [path] output table file\n";
	print "\t-1 [switch] skip first column\n";
	print "\t-2 [switch] skip second column\n";
	print "\t-l [file:name_col#:len_col#] inputs length file for speed (default columns 0:1";
	print " It requires standard miropeat output name b e len name2 b2 e2 len2\n";
	exit "***********************************************************\n";
}

getopts('t:o:s:f:cl:12');
$opt_t || die "Please enter with -t the path for the pairwise blast table\n";
if (!$opt_c) {
	$opt_f || die "Please enter with -f the path to the unfractionated files\n";
}
$opt_s || die "Please enter with -s the size of the fragments\n";
$opt_o ||= $opt_t ."\.defrac";



if ($opt_l) {
	my ($file,$col1,$col2)=split "\t", $opt_l;
	$col1=0 if !defined $col1;
	$col2=1 if !defined $col2;
	open ( INLIST, "$file") || die "Can't read input length list ($file)\n";
	while (<INLIST>) {
		next if !/\t/;
		my @c= split "\t";
		$len{$c[$col1]}=$c[$col2];
	}
	close INLIST;
	

}



open (TABLE, $opt_t) || die "Can't open pairwise blast table: $opt_t\n";
open (NEWTABLE, ">$opt_o") || die "Can't output table file: $opt_o\n";

my $table_header=<TABLE>;
print NEWTABLE "$table_header";
$count=0;
while (<TABLE>) {
	s/\r\n/\n/;
	chomp;
	my @c=split "\t";
	#####GET ORIGINAL LENGTH#####
		my $name=$c[0];
	if ($c[0]=~/(.*)_(\d+)$/ && !$opt_1) {
		my ($name,$part)=($1,$2);
		$c[0]=$name;
		$c[1]+=$part*$opt_s;
		$c[2]+=$part*$opt_s;
		$len{$name}='NA' if $opt_f eq 'SKIP';
		if (!defined $len{$name} && !$opt_c) {
			my $filepath=&find_file ( $opt_f,"",$name,"");
			$len{$name}=&get_length("$filepath");
		}
		$c[3]=$len{$name} if !$opt_c;
	}
	if ($c[4]=~/(.*)_(\d+)$/ && !$opt_2) {
		my ($name,$part)=($1,$2);
		$c[4]=$name;
		$c[5]+=$part*$opt_s;
		$c[6]+=$part*$opt_s;	
		$len{$name}='NA' if $opt_f eq 'SKIP';
		if (!defined $len{$name} & !$opt_c) {
			my $filepath=&find_file ( $opt_f,"",$name,"");
			$len{$name}=&get_length("$filepath");
		}
		$c[7]=$len{$name} if !$opt_c;
	}
	
	print NEWTABLE join("\t",@c)  , "\n";
}


sub find_file {
	my ($paths, $sub_paths, $names, $extensions) = @_;
	my @paths=split ":",$paths;
	my @sub_paths=split ":",$sub_paths;
	@sub_paths=('') if !@sub_paths;
	my @names=split ":",$names;
	my @extensions=split ":",$extensions;
	@extensions=('') if ! @extensions;
	for my $path (@paths) {
	 for my $sub_path (@sub_paths) {
	 	for my $name (@names) {
	 		for my $ext (@extensions) {	
	 			my $p = $path;
	 			$p .= "/$sub_path" if $sub_path;
	 			$p .= "/$name$ext";
	 			print "FILE $p exists? ";
	 			if (-e $p) {print "YES\n"; return $p; }
	 			print "NO\n";
	 		}
	 	}
	 }
	}
	print "    COULD NOT FIND FILE\n";
	return "";
}

sub get_length {
		open (FASTA,"$_[0]") || die "Can not open fasta file $_[0]\n";
		my $header=<FASTA>;
		my ($totbases,$a,$g,$c,$t,$n);
		die ("File does not contain a fasta header!\n") if $header !~/^>/;
		while ( <FASTA> ) {
			s/\r\n/\n/;
			s/ +//;
			chomp;
			chomp;
			$_ = uc($_);
			$totbases += length;
			$a+= tr/A/A/;
			$g+= tr/G/G/;
			$c+= tr/C/C/;
			$t+= tr/T/T/;
			$n+= tr/N/N/;
	
		}
		print "TOTAL BASE PAIRS: $totbases A:$a C:$c G:$g T:$t N:$n O:",$totbases-$a-$c-$g-$t-$n,"\n";
		return $totbases;
}