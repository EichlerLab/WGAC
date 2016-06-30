#!/usr/bin/perl

$tableFile=$ARGV[0];
$fastaDir =$ARGV[1]; 
$outputDir=$ARGV[2];

$fileName=$tableFile;
$fileName=~ s/^.*\///;

$workDir ="/tmp/endtrim";

$progDir ="/net/eichler/vol5/home/ssajjadi/wgacbin/step_8_mpi";

open(OUT,">>$outputDir/output.txt");

system "perl $progDir/pairwise.pl -x 0.3:3000:-3:-20 -n -f $fastaDir -t $tableFile  -u $fileName:random"; 
print OUT "perl $progDir/pairwise.pl -x 0.3:3000:-3:-20 -n -f $fastaDir -t $tableFile  -u $fileName:random\n";

$tmpDir=$fileName;
$tmpDir.="_tmp";
#system "rm -rf $workDir/$tmpDir/*";

system "perl $progDir/pairwise.pl -x 0.3:3000:-3:-20 -c -f $fastaDir -t $tableFile -u $fileName:random"; #>output: output8 (tableFile.trim)
print OUT "perl $progDir/pairwise.pl -x 0.3:3000:-3:-20 -c -f $fastaDir -t $tableFile -u $fileName:random\n"; #>output: output8 (tableFile.trim)

system "mkdir $outputDir/$fileName";
system "cp -rf $workDir/$fileName.trim $outputDir/$fileName";
system "rm -rf $workDir/$tmpDir";
system "rm -rf $workDir/localTrim_$tmpDir";
system "rm -rf $workDir/$fileName.trim";
