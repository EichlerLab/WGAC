

#!/usr/bin/perl

## this is used to filter out the redundancy in genomicSuperDup
## after the superdup is generated

use strict;
my @a;
my %hash;
my $id;
my $value;

while(<>){
		@a = split(/\s+/);
		$id = "$a[0]$a[1]$a[2]";
		$value = "$a[3]$a[4]$a[5]";

#	 print "$id and -------  $value\n";
	
	if ($hash{$value}eq $id){
#			print " find ------------ has $id-$hash{$id}-$value\n";
	}else{
			print $_;
			$hash{$id} = $value;
	}

}
	
			
