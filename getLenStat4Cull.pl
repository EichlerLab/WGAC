

#!/usr/bin/perl
use strict;
my @a;
my $intraCount = 0;
my $interCount = 0;
my ($intraLenq, $interLenq, $intraLens, $interLens, $plen);


## first to set up variables
my %interLenCountHash = (
	1 => 0,
	2 => 0,
	3 => 0,
	4 => 0,
	5 => 0,
	6 => 0,
	7 => 0,
	8 => 0,
	9 => 0,
	10 => 0,
	20 => 0,
	30 => 0,
	40 => 0,
	50 => 0,
);

my %intraLenCountHash = %interLenCountHash;

my %interLenHash = (
  1 => 0,
  2 => 0,
  3 => 0,
  4 => 0,
  5 => 0,
  6 => 0,
  7 => 0,
  8 => 0,
  9 => 0,
  10 => 0,
  20 => 0,
  30 => 0,
  40 => 0,
  50 => 0,
);

my %intraLenHash = %interLenHash;
my ($interLenS, $intraLenS);
while(<>){
  if(/QNAME/){
#			print $_;
			next;
	}    
	if (/_random/){
		s/_random//g;
	}

	@a = split (/\t/);
	$plen = $a[22];
	if($a[0] eq $a[4]){		# intra chrom dup
			$intraCount++;
			$intraLenq = abs($a[1]-$a[2]);
			$intraLenS = abs($a[5]-$a[6]);
			sortLength($intraLenq, "intra");
			sortLength($intraLenS, "intra");
	}else{								#inter
			$interCount++;
			$interLenq = abs($a[1]-$a[2]);
			$interLenS = abs($a[5]-$a[6]);
			sortLength($interLenq,"inter");
			sortLength($interLenS,"inter");
	}
}
print "inter pairs is $interCount\n";
print "intra pairs is $intraCount\n";
##-- to print out the intra length distribution and interlength distributions
my @keyArray = keys(%interLenCountHash);
@keyArray = sort {$a <=> $b } @keyArray;
my $k;
my $total = 0;
my $tlen =0;
print "size\tinter\tintra\tTotal\tinterlen\tintralen\ttlen\n";
foreach $k(@keyArray){
	$total = $interLenCountHash{$k} + $intraLenCountHash{$k};
	$tlen = $interLenHash{$k} + $intraLenHash{$k};
	print "$k.kb\t$interLenCountHash{$k}\t$intraLenCountHash{$k}\t$total\t$interLenHash{$k}\t$intraLenHash{$k}\t$tlen\n";
}
#-- subs ---- 
sub sortLength{
	my($len, $type) = @_;
	if ($type eq "inter"){
		if ($len > 999 && $len <2000){
				$interLenCountHash{1}++;
				$interLenHash{1} += $len;
		}elsif ($len > 1999 && $len < 3000){
				$interLenCountHash{2}++;
				$interLenHash{2} += $len;
		}elsif ($len > 2999 && $len < 4000){
        $interLenCountHash{3}++;
				$interLenHash{3} += $len;	
		}elsif ($len > 3999 && $len < 5000){
        $interLenCountHash{4}++;
			$interLenHash{4} += $len;
		}elsif ($len > 4999 && $len < 6000){
        $interLenCountHash{5}++;
				$interLenHash{5} += $len;
		}elsif ($len > 5999 && $len < 7000){
        $interLenCountHash{6}++;
				$interLenHash{6} += $len;
		}elsif ($len > 6999 && $len < 8000){
        $interLenCountHash{7}++;
				$interLenHash{7} += $len;
		}elsif ($len > 7999 && $len < 9000){
        $interLenCountHash{8}++;
				$interLenHash{8} += $len;
		}elsif ($len > 8999 && $len < 10000){
        $interLenCountHash{9}++;
				$interLenHash{9} += $len;
		}elsif ($len > 9999 && $len < 20000){
        $interLenCountHash{10}++;
				$interLenHash{10} += $len;
		}elsif ($len > 19999 && $len < 30000){
        $interLenCountHash{20}++;
				$interLenHash{20} += $len;
		}elsif ($len > 29999 && $len < 40000){
        $interLenCountHash{30}++;
				$interLenHash{30} += $len;
		}elsif ($len > 39999 && $len < 50000){
        $interLenCountHash{40}++;
				$interLenHash{40} += $len;
		}elsif ($len > 49999 && $len < 4000000){
        $interLenCountHash{50}++;
				$interLenHash{50} += $len;
		}

	}elsif($type eq "intra"){
		if ($len > 999 && $len <2000){
        $intraLenCountHash{1}++;
				$intraLenHash{1} += $len;
    }elsif ($len > 1999 && $len < 3000){
        $intraLenCountHash{2}++;
				$intraLenHash{2} += $len;
    }elsif ($len > 2999 && $len < 4000){
        $intraLenCountHash{3}++;
				$intraLenHash{3} += $len;
    }elsif ($len > 3999 && $len < 5000){
        $intraLenCountHash{4}++;
				$intraLenHash{4} += $len;
    }elsif ($len > 4999 && $len < 6000){
        $intraLenCountHash{5}++;
				$intraLenHash{5} += $len;
    }elsif ($len > 5999 && $len < 7000){
        $intraLenCountHash{6}++;
				$intraLenHash{6} += $len;
    }elsif ($len > 6999 && $len < 8000){
        $intraLenCountHash{7}++;
				$intraLenHash{7} += $len;
    }elsif ($len > 7999 && $len < 9000){
        $intraLenCountHash{8}++;
				$intraLenHash{8} += $len;
    }elsif ($len > 8999 && $len < 10000){
        $intraLenCountHash{9}++;
				$intraLenHash{9} += $len;
    }elsif ($len > 9999 && $len < 20000){
        $intraLenCountHash{10}++;
				$intraLenHash{10} += $len;
    }elsif ($len > 19999 && $len < 30000){
        $intraLenCountHash{20}++;
				$intraLenHash{20} += $len;
    }elsif ($len > 29999 && $len < 40000){
        $intraLenCountHash{30}++;
				$intraLenHash{30} += $len;
    }elsif ($len > 39999 && $len < 50000){
        $intraLenCountHash{40}++;
				$intraLenHash{40} += $len;
    }elsif ($len > 49999 && $len < 4000000){
        $intraLenCountHash{50}++;
				$intraLenHash{50} += $len;
    }

	}
}
