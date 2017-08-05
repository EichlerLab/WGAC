

#!/usr/bin/perl
use strict;
my @a;
my $intraCount = 0;
my $interCount = 0;
my ($intraSim, $interSim, $plen);


## first to set up variables
my %interSimCountHash = (
	0.90 => 0,
	0.91 => 0,
	0.92 => 0,
	0.93 => 0,
	0.94 => 0,
	0.95 => 0,
	0.96 => 0,
	0.97 => 0,
	0.98 => 0,
	0.99 => 0,
	0.995 => 0,
	1 => 0
);

my %intraSimCountHash = %interSimCountHash;

my %interSimLengthHash = (
  0.90 => 0,
  0.91 => 0,
  0.92 => 0,
  0.93 => 0,
  0.94 => 0,
  0.95 => 0,
  0.96 => 0,
  0.97 => 0,
  0.98 => 0,
  0.99 => 0,
  0.995 => 0,
  1 => 0
);

my %intraSimLengthHash = %interSimLengthHash;

while(<>){
  if(/QNAME/){
#			print $_;
			next;
	}    
	if(/random/){
		s/_random//g;
	}
	
	@a = split (/\t/);
	$plen = $a[22];
	if($a[0] eq $a[4]){		# intra chrom dup
			$intraCount++;
			$intraSim = $a[27];
			sortSim($intraSim, "intra", $plen);
			sortSim($intraSim, "intra", $plen);
	}else{								#inter
			$interCount++;
			$interSim =$a[27];
			sortSim($interSim,"inter", $plen);
			sortSim($interSim,"inter", $plen);
	}
}
print "inter pairs is $interCount\n";
print "intra pairs is $intraCount\n";
##-- to print out the intra length distribution and interlength distributions
my @keyArray = keys(%interSimCountHash);
@keyArray = sort {$a <=> $b } @keyArray;
my $k;
my $total = 0;
my $tlen =0;
print "size\tinter\tintra\tTotal\interLength\tintraLength\tTotalLength\n";
foreach $k(@keyArray){
	$total = $interSimCountHash{$k} + $intraSimCountHash{$k};
	$tlen = $interSimLengthHash{$k} + $intraSimLengthHash{$k};
	print "$k\t$interSimCountHash{$k}\t$intraSimCountHash{$k}\t$total\t$interSimLengthHash{$k}\t$intraSimLengthHash{$k}\t$tlen\n";
}
#-- subs ---- 
sub sortSim{
	my($sim, $type, $len) = @_;
	if ($type eq "inter"){
		if ($sim >= 0.9 && $sim < 0.91){
				$interSimCountHash{0.9}++;
				$interSimLengthHash{0.9} += $len;
		}elsif ($sim >= 0.91 && $sim < 0.92){
				$interSimCountHash{0.91}++;
				$interSimLengthHash{0.91} += $len;
		}elsif ($sim >= 0.92 && $sim < 0.93){
        $interSimCountHash{0.92}++;
				$interSimLengthHash{0.92} += $len;
		}elsif ($sim > 0.93 && $sim < 0.94){
        $interSimCountHash{0.93}++;
				$interSimLengthHash{0.93} += $len;
		}elsif ($sim > 0.94 && $sim < 0.95){
        $interSimCountHash{0.94}++;
				$interSimLengthHash{0.94} += $len;
		}elsif ($sim > 0.95 && $sim < 0.96){
        $interSimCountHash{0.95}++;
				$interSimLengthHash{0.95} += $len;
		}elsif ($sim > 0.96 && $sim < 0.97){
        $interSimCountHash{0.96}++;
				$interSimLengthHash{0.96} += $len;
		}elsif ($sim > 0.97 && $sim < 0.98){
        $interSimCountHash{0.97}++;
				$interSimLengthHash{0.97} += $len;
		}elsif ($sim > 0.98 && $sim < 0.99){
        $interSimCountHash{0.98}++;
				$interSimLengthHash{0.98} += $len;
		}elsif ($sim > 0.99 && $sim < 0.995){
        $interSimCountHash{0.99}++;
				$interSimLengthHash{0.99} += $len;
		}elsif ($sim > 0.995 && $sim <= 1){
				if ($sim == 1){
				$interSimCountHash{1}++;
				$interSimLengthHash{1} += $len;
				}else{
        $interSimCountHash{0.995}++;
				$interSimLengthHash{0.995} += $len;
				}
		}
	}elsif($type eq "intra"){
	 if ($sim >= 0.9 && $sim < 0.91){
        $intraSimCountHash{0.9}++;
				$intraSimLengthHash{0.9} += $len;
    }elsif ($sim >= 0.91 && $sim < 0.92){
        $intraSimCountHash{0.91}++;
				$intraSimLengthHash{0.91} += $len;
    }elsif ($sim >= 0.92 && $sim < 0.93){
        $intraSimCountHash{0.92}++;
				$intraSimLengthHash{0.92} += $len;
    }elsif ($sim > 0.93 && $sim < 0.94){
        $intraSimCountHash{0.93}++;
				$intraSimLengthHash{0.93} += $len;
    }elsif ($sim > 0.94 && $sim < 0.95){
        $intraSimCountHash{0.94}++;
				$intraSimLengthHash{0.94} += $len;
    }elsif ($sim > 0.95 && $sim < 0.96){
        $intraSimCountHash{0.95}++;
				$intraSimLengthHash{0.95} += $len;
    }elsif ($sim > 0.96 && $sim < 0.97){
        $intraSimCountHash{0.96}++;
				$intraSimLengthHash{0.96} += $len;
    }elsif ($sim > 0.97 && $sim < 0.98){
        $intraSimCountHash{0.97}++;
				$intraSimLengthHash{0.97} += $len;
    }elsif ($sim > 0.98 && $sim < 0.99){
        $intraSimCountHash{0.98}++;
				$intraSimLengthHash{0.98} += $len;
    }elsif ($sim > 0.99 && $sim < 0.995){
        $intraSimCountHash{0.99}++;
				$intraSimLengthHash{0.99} += $len;
    }elsif ($sim > 0.995 && $sim <= 1){
        if ($sim == 1){
        $intraSimCountHash{1}++;
				$intraSimLengthHash{1} += $len;
        }else{
        $intraSimCountHash{0.995}++;
				$intraSimLengthHash{0.995} += $len;
        }
    }


	}
}
