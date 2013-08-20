#!/usr/bin/perl

## this is used to generate the mask out file from a masked fasta file
## based on the lower case masking

use strict;
my @a;
my $uCount = 0;
my $chrCount = 0;
my $tigCount = 0;
my $def ;
my $line;
my $tigbp = 0;
my $ubp = 0;
my $outfile = $ARGV[1];
my $infile = $ARGV[0];

open (O, ">$outfile");
## and print the header, t=for 2 line and a space
print O " SW  perc perc perc  query  position in query  matching  repeat  position in  repe at\n";
print O "score  div. del. ins.  sequence  begin  end (left)  repeat  class/family  begin  end(left)ID\n\n";

open(F, "$infile");
while(<F>){
      chomp;
	if (/>/){
		s/\>//;
		$def = $_;
	}else{
		$line .= $_;
		$uCount++;
	}
}
processNprint(\$line);

close(F);
close(O);

#### subs ---------------------

sub processNprint{
	my $refline = shift;
	my @seq = split(//, $$refline); # individual letters
	my $len = @seq;
	my $s;
	my $pos = 0;
	my $lowCaseFlag = 0;
	my $lowCaseStart = 0;
	my $lowCaseEnd = 0;
	my $order = 0;

	foreach $s(@seq) {
		$pos++; # 1-based coordinates
#		print "$s ---- \n";
		if ($s =~/[N]/){
				# Ns
				if($lowCaseFlag == 0){	# start
					$lowCaseFlag = 1;
					$lowCaseStart = $pos;
					$order++;
				}

		} else { # else it's uppercase
				if ($lowCaseFlag == 1){		# the end of the lowcase
					$lowCaseFlag = 0;
					$lowCaseEnd = $pos - 1;
					print O "1000  20 3.9  2.0  $def  $lowCaseStart  $lowCaseEnd ($len) +  repeat  class/Fam    1  2000 (100) $order\n";
				}

		}
	}

	# if sequence ends in lowercase letter
	if ($lowCaseFlag == 1) {
		$lowCaseEnd = $pos;
		print O "1000  20 3.9  2.0  $def  $lowCaseStart  $lowCaseEnd ($len) +  repeat  class/Fam  1  2000 (100) $order\n";
	}

}
