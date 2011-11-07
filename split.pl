

#!/usr/bin/perl
## thisis used for spliting the file 
## and create a dir
# and for each file mbeing put into a dir


use strict;
my @a;
my $file = $ARGV[0];					# input files
my $numberToSplit = $ARGV[1];	#number of dir to split to
my $outfileName = 0;
my $count = 0;
my $fileCount = 1;
my $d = "newdir";
my $header = "header";
my $tLine= $ARGV[2];

if(!$ARGV[2]){
	print " this is use to split the both.parse.defugu file\n";
	print "USAGE: script both.parse.defugu numberofFiles totalLineinBoth.parse.defugu -- \n";
	die;
}

print "$tLine\n";


my $linePerFile = $tLine/$numberToSplit + 2;  # to make sure to get all lines
print "the line per file is $linePerFile\n";
print "total count is $tLine\n";
print "$numberToSplit---\n";
## to create a dir
system "mkdir newdir";
open (F, "$file");
while(<F>){
  if(/QNAME/){			# header
				$d = "newdir/$fileCount";
				system "mkdir $d";
				$outfileName = $fileCount;
			 	open (O, ">$d/$outfileName");
        print O "$header";
        print O "$_";
        $count = 1;

			$header = $_;
	}else{
			if ($count > $linePerFile){
				close(O);
				$fileCount++;
				$d = "newdir/$fileCount";
        system "mkdir $d";
        $outfileName = $fileCount;
				print "generate $fileCount files and  line per file is $linePerFile\n";
				open (O, ">$d/$outfileName");
				print O "$header";
				print O "$_";
				$count = 1;
			}else{
					print O "$_";
					$count++;
			}

	}
}


			
