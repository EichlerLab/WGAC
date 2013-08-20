#!/usr/bin/env perl
#10/01/12 Tin Louie changed rm -rf to rmdir (allow another processor to finish using the subdirectory)
#10/01/08 Tin Louie added test of the exit code of mkdir (test if another processor is working on current line)
#10/01/06 Tin Louie modified calls to find_file (to look for files with 4-digit extensions)
#01/02/30 adding ability to 
#01/01/09 added program stats#
#00/09/22 allow multiple output directories###
#00/05/29 Major modifications to make it batch parrallel #####
#      add
#00/13/14 4 PM mofified to take multiple directory input and find the file
#      added find file subroutine and -f option multiples with :
#		 modified &find_file for empty extensions and empty subdirectories
#		 added -r option for starting row number when not using -u

use strict 'vars';
use Getopt::Std;
use FindBin;

use vars qw($program $pversion $pdescription $pgenerate $pusage);
$program = $0;
$program =~ s/^.*\///;
### program stats ###
$pversion='3.010109';
$pdescription = "$program (ver:$pversion) generates global alignments using align_fast2.pl on a table of pairwise coordinates allowing for multiple processors";
$pgenerate= 'total: labmates: genetics: public:';
$pusage="$program -t [table] -i [col f:b:e] -j [col f:b:e] -f [fasta directories] ..";
### program stats end ###

use vars qw($opt_t $opt_i $opt_j $opt_f $opt_x $opt_d $opt_p $opt_u $opt_b $opt_o $opt_n $opt_s $opt_l $opt_r);
use vars qw(@header $temp $row $ppath $tmpdir $puid $odnumber);

if (!defined $ARGV[0]) { 
	print "USAGE\n";
	print "$pusage\n";
	print "DESCRIPTION\n$pdescription\n";
	print "REQUIRED ARGUMENTS\n";
	print "-t  input table from which to draw data\n";
	print "-i  input sequence1 columns [fasta1:begin1:end1]\n";
	print "-j  input sequence2 columns [fasta2:begin2:end2]\n";
	print "-f  [dir path(s):]  directories containing fasta files\n";
	print "      (separate multiple paths with colon(:)\n";
	print "-u  [buid:puid] batch unique id and processor unique id : separated\n";
	print "-o  [path] directory to storeout files\n";
	print "OPTIONAL ARGUMENTS:\n";
	print "-l  [bp] size of pieces to align (default 60000)\n";
	print "-d  [integer] number of alignments per directory\n";
	print "-s  [column] sufix unique id for each alignment (default use row #)\n";
	print "-n  [column] prefix name for alignments\n";
	print "-x  [column] containing the name for the file\n";
	print "-b  [integer] number of bases for fractionation of fasta files (default no fractionatition)!\n";
	print "-r  [integer] starting row number for alignment naming (default=1)\n";
	print "-p  [bp] add additional sequence to ends for robust overlap \n";

	exit ;
}
$ppath="";

$ppath="perl /binp/" if $^O eq "MSWin32";
###GET OPTIONS
getopts('t:i:j:p:f:o:n:j:u:l:r:s:d:x:b:');
print STDERR "INPUT FASTAS AT $opt_f\n";
open (TABLE, "$opt_t") || die "Can't open table $opt_t\n!";
$temp = <TABLE>; #header row#
$opt_l ||= 60000;
$opt_r ||=1;
$opt_b ||='';
if ($opt_b) {
	$opt_b = "-b $opt_b";
	print STDERR "$opt_b\n";
}
$opt_d ||= 9999999999999;
$opt_u ||='alignfast:random';
$opt_x ||= '';
$opt_o || die "Please input with -o a directory to store alignments in!\n";
$odnumber='0000';

#SETUP BUID AND PUID#
die "Format of -u must be buid:puid\n" if $opt_u !~ /\:/;
($tmpdir,$puid) = split ":", $opt_u;
if ($puid eq 'random') {
 $puid = "S" . int (rand(100000000000));
}
$tmpdir.="_tmp";
system "mkdir $tmpdir";
system "mkdir $opt_o";
system "mkdir $opt_o/$odnumber";
my $maindir=`pwd`;
chomp $maindir;
$maindir=~s/\/$//;
####MAIN LOOP ######
$row=$opt_r-1;
while (<TABLE>) {
	$row++;
	my $newdirnum=int (($row-0.8)/$opt_d);
	print STDERR "ROW:$row  NEWDIRNUM $newdirnum  ODNUM:$odnumber\n";
	if ($newdirnum > $odnumber) {
		$odnumber=substr ("0000".$newdirnum,-4);
		system "mkdir $opt_o/$odnumber";
		print STDERR "MAKING NEW DIRECTORY $opt_o/$odnumber\n";
	}
	print STDERR "######ROW $row ... \n";
	chomp;
	my @col = split "\t";
	
	#DEFINE NAME#
	my $name=$opt_n;
	#Changed by Saba on 26April 2009 from substr("000000".$row,-6) to substr("0000000".$row,-7)
	if (defined $opt_s) {$name.=$col[$opt_s];} else {$name.=substr("0000000".$row,-7);}
	if ($opt_x ne '') {
		$name=$col[$opt_x];
	}

	##TEST FOR EXISTENCE#
	if (-e "$opt_o/$odnumber/$name") {
		print STDERR "    ...already exists in ($opt_o/$odnumber).\n";
		next;
	}
	my $exitcode = system "mkdir $tmpdir/$name";
	if ($exitcode != 0) {
		print STDERR "    ...another processor is already working on it.\n";
		next;
	}
	
	##EXTRACT COLUMN INFORMATION FROM TABLE####
	my ($s1,$b1,$e1)=('','','');
	 if ($opt_i =~ /^(\d+)/ ) {
		$s1=$1;
		($b1,$e1)=($1,$2) if $opt_i=~ /^\d+:(\d+):(\d+)$/;
	} else {
			 	die "-i ($opt_i) is improper format column# seq:begin(opt):end(opt)\n";
	}
	$s1=$col[$s1] if  $s1 ne '';
	$b1=$col[$b1] if  $b1 ne '';
	$e1=$col[$e1] if  $e1 ne '';
	my ($s2,$b2,$e2)=('','','');
	if ($opt_j =~/^(\d+)/ ) {
		$s2=$1;
		($b2,$e2)=($1,$2) if $opt_j=~ /^\d+:(\d+):(\d+)$/;
	} else {
	 	die "-j ($opt_j) is improper format column# seq:begin(opt):end(opt)\n";
	
	}
	$s2=$col[$s2] if $s2 ne '';
	$b2=$col[$b2] if $b2 ne '';
	$e2=$col[$e2] if $e2 ne '';
	print STDERR "$s1 $b1-$e1  $s2 $b2-$e2\n";
	##CHECK FOR FASTA FILES#
	my $path1 = &find_file($opt_f, '',$s1, ':_000:_0000');
	if ($opt_b) {
		$path1 =~ s/_0{3,4}$//;  # remove 3 or 4 zeroes
	}
	print STDERR "FOUND $path1\n";
	if (!$path1 ) {
		print STDERR "S1:$s1 doesn't exist--Skipping ($opt_f $s1)\n";
		next;
	}
	
	my $path2 = &find_file($opt_f, '',$s2, ':_000:_0000');
	if ($opt_b) {
		$path2 =~ s/_0{3,4}$//; # remove 3 or 4 zeroes
	}
	print STDERR "FOUND $path2\n";
	if (!$path2) {
		print STDERR "S2:$s2 doesn't exist--Skipping ($opt_f)\n";
		next;
	}

	my ($frag1, $frag2);
	if ($b1 eq '' && $b2 eq '' && $e2 eq '' && $e1 eq '') {
		$frag1='';
		$frag2='';
	} else {
	
		if ($b1 eq "" || $b2 eq "" || $e1 eq "" || $e2 eq "") {
			print STDERR "Missing numbers for b1 $b1 b2 $b2 e1 $e1 e2 $e2\n";
			next;
		}
		###ADD EXTRA BASES TO AN ALIGNMENT###
		if ($opt_p ne "") {
			my ($plus,$l1,$l2) = ($1,$col[$2],$col[$3]);
			$b1-= $plus;
			$b1 =1 if $b1 < 1;
			$e1+= $plus;
			$e1=$l1 if $e1 > $l1;
			print "L1 $l1   L2 $l2\n";
			if ($b2 < $e2) {
				$b2-= $plus;
				$b2 =1 if $b2 < 1;
				$e2+= $plus;
				$e2 = $l2 if $e2 > $l2;
			} else {
				$b2+= $plus;
				$b2 = $l2 if $b2 > $l2;
				$e2-= $plus;
				$e2 =1 if $e2 < 1;
			}
		}
		$frag1=":$b1:$e1";
		$frag2=":$b2:$e2";
	}
	
	print STDERR "#RUN ALIGN_FASTA3.PL#\n";
	#chdir "$tmpdir/$name";
	#system "ls";
	my $ofilename="$maindir/$opt_o/$odnumber/$name";
	my $command= "$FindBin::Bin/align_fast3.pl -i $path1$frag1 -j $path2$frag2 -l $opt_l -o $ofilename -f -40 -g -1 $opt_b";
	#eray changed removed the maindirs before path1 and path2
	print "mkdir -p $maindir/$tmpdir/$name;cd $maindir/$tmpdir/$name;$command;cd $maindir;rm -rf $maindir/$tmpdir/$name\n";  
	print STDERR "COMMAND $command\n";
	

	

	# allow another processor to finish using tmpdir/name, if files exist within this subdirectory
	#system "rm -rf $tmpdir/$name";
	#system "rmdir --ignore-fail-on-non-empty $tmpdir/$name";
	 
	
}


sub find_file {
	my ($paths, $sub_paths, $names, $extensions) = @_;
	
	my @paths=split ":",$paths;
	my @sub_paths=split ":",$sub_paths;
	@sub_paths=("") if !@sub_paths;
	my @names=split ":",$names;
	my @extensions=split /:/,$extensions;
	@extensions=("") if !@extensions;
	print STDERR join ("XXX", @extensions), "\n";
	for my $path (@paths) {
	 for my $sub_path (@sub_paths) {
	 	for my $name (@names) {
	 		for my $ext (@extensions) {	
	 			my $p = $path;
	 			$p .= "/$sub_path" if $sub_path;
	 			$p .= "/$name$ext";
	 			#my $wpd=`wpd`;
	 			print STDERR "TESTING $p\n";
	 		
	 			return $p if (-e $p);
	 		}
	 	}
	 }
	}
	my $pause=<STDIN>;
	return "";
}




