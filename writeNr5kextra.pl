

#!/usr/bin/perl


## this is used to generate extra files for Orangutan dup overlay
## on human chromsome
## for extra file generation
## input files are those splitted dup cordinants
## output into a new directory
## make sure to check the files and change directory name after finished
##
## $0 dir
##
## lin chen 2007-3-6

use strict;

my @a;
my $fileNames = $ARGV[0];
my $a;
my $dif;
print "seqName\tbegin\tend\tcolor\toffset\twidth\n";
while(<>){
      chomp;
			@a = split/\t/;
		$dif = abs($a[1] - $a[2]);
			if ($dif >= 5000){
#			print "$_\n";
#			print "$_\tbrown\n";
			print "$a[0]\t$a[1]\t$a[2]\tblack\t8\t6\t\twgac\n";
			}
}			

