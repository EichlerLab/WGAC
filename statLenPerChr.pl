

#!/usr/bin/perl

# this is use to count the length in Merge file
#-- list for each chromosome
use strict;

my @a;
my %chrLen;
my %hitLen;
my $a = "chr";
my $total = 0;

while(<>){
      
		chomp;
		@a = split (/\s+/);
		if ($hitLen{$a[0]}){
			$hitLen{$a[0]} = $hitLen{$a[0]} + $a[3];
		
		}else{
				$hitLen{$a[0]} = $a[3];
		}
		$chrLen{$a[0]} = $a[4];
	

}
@a = keys(%hitLen);

@a = sort @a;
my $percent;
foreach $a(@a){
	if($chrLen{$a} == 0){
			print "$a and $hitLen{$a}===\n";
	}else{
		$percent = $hitLen{$a}/$chrLen{$a};
		$total = $total + $hitLen{$a};
	#	print "$a -----\n";
		print "$a\t$hitLen{$a}\t$chrLen{$a}\t$percent\n";
	}
}

print "Total\t$total\n";
			
