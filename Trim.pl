#!/usr/bin/perl


$tableFile=$ARGV[0];
$fastaDir =$ARGV[1];
$outputDir=$ARGV[2];
$species=$ARGV[3];
$usrname=`whoami`;
chomp $usrname;
$fileName=$tableFile;
$fileName=~ s/^.*\///;
$workDir="/tmp/$usrname/$species/endtrim";

# changed from /net/eichler/vol4/home/jlhudd/wgac by DG (8/16/2019)
$progDir = `pwd`;
chomp( $progDir );

print "outputDir = $outputDir\n";

open(OUT,">>$outputDir/output.txt");


system "perl $progDir/pairwise.pl -t $tableFile -f $fastaDir -u $fileName:random -d $workDir -x 0.3:3000:-3:-20 -n";
print OUT "perl $progDir/pairwise.pl -t $tableFile -f $fastaDir -u $fileName:random -d $workDir -x 0.3:3000:-3:-20 -n\n";

$tmpDir=$fileName;
$tmpDir.="_tmp";


system "perl $progDir/pairwise.pl -t $tableFile -f $fastaDir -c -u $fileName:random -d $workDir -x 0.3:3000:-3:-20"; #>output: output8 (tableFile.trim)
print OUT "just executed: perl $progDir/pairwise.pl -t $tableFile -f $fastaDir -c -u $fileName:random -d $workDir -x 0.3:3000:-3:-20\n"; #>output: output8 (tableFile.trim)

system "mkdir -p $outputDir/$fileName";
system "cp -rf $workDir/$fileName.trim $outputDir/$fileName";

system "rm -rf $workDir/$tmpDir";
system "rm -rf $workDir/localTrim_$tmpDir";
system "rm -rf $workDir/$fileName.trim";
