

#!/usr/bin/perl


## this is use to perform batch work for a directory
## to convert all the ps to pdf file
## $0 dir
##
## lin chen 2007-3-6

use strict;

my (@a, @files);
my $d = $ARGV[0];	# the directory
my $fileName = "file name";
my $dirNames = "sub dir name";
my $a;
my $f;
my $outfileName;


$dirNames = `ls $d`;
# print "$fileNames\n";
@a =split(/\n/, $dirNames);
foreach $a(@a){
	if($a =~ /starb/){
	
		$fileName = `ls $d/$a`;
		@files = split(/\n/, $fileName);
		
		foreach $f (@files){
			$outfileName = $f;
			$outfileName =~ s/\.ps/\.pdf/;
			print " doing ps2pdf $d/$a/$f $d/$a/$outfileName -- \n";
			system "ps2pdf $d/$a/$f $d/$a/$outfileName";
			print "processed from $d/$a/$f to $d/$a/$outfileName\n";
		}
	}
	
}	

