#!/usr/bin/perl
## thisis used for combining the culls in the last weild step 
## to make a single join.cull table
## since the process need to have file splitted 
## after the process , i need to put them back togather


use strict;
my @a;
my $d = $ARGV[0];

my $outputfile ="$d/ParallelOutput.trim"; #"$d/contig.join.all.cull";
my $inputfile = "file path and name";

my $count = 0;

open(O, ">$outputfile");
opendir ( DIR, $d ) || die "Error in opening dir $d\n";

my @subDirList = readdir(DIR);

#my $it;

#foreach $it(0..$#subDirList){
#	$subDirList[$it] =~ s/^(part)(\d+)$/\2/;
#}
my @subDirList = sort {$a <=> $b} @subDirList;

my $i = 1;
my $subDir;
foreach $subDir(@subDirList){
	if($subDir=~/[^0-9]/){
		next;
	}
	#$subDir="part$subDir";
	$inputfile = "$d/$subDir/$subDir.trim";#"$d/$subDir/oo.weld10kb.$subDir.join.all.cull";
	print "the input file $inputfile\n";
	#print "the output file is $outputfile\n";
	open(F, "$inputfile")|| die "can not open the file $inputfile\n";;
	while (<F>){
		if (/QNAME/){
			if ($i == 1){
				s/headerQNAME/QNAME/;
				print O "$_";				
			}
		}else{		
			print O "$_";
			$count++
		}
	}
	$i++;
	close(F);
	#print "$count lines has been written to $outputfile\n";
	$count = 0;
}
close (O);
