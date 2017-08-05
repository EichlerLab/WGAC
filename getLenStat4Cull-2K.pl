

#!/usr/bin/perl
use strict;
my @a;
my $intraCount = 0;
my $interCount = 0;
my ($intraLenq, $interLenq, $intraLens, $interLens, $plen);

my %nameHash = (
	1 => "1-2",
	2 => "2-4",
	4 => "4-6",
	6 => "6-8",
	8 => "8-10",
	10 => "10-15",
	15 => "15-20",
	20 => "20-25",
	25 => "25-30",
	30 => "30-35",
	35 => "35-40",
	40 => "40-45",
	45 => "45-50",
	50 => "50-60",
	60 => "60-70",
	70 => "70-80",
	80 => "80-90",
	90 => "90-100",
	100 => ">100",

);

## first to set up variables
my %interLenCountHash = (
	1 => 0,
	2 => 0,
	4 => 0,
	6 => 0,
	8 => 0,
	10 => 0,
	15 => 0,
	20 => 0,
	25 => 0,
	30 => 0,
	35 => 0,
	40 => 0,
	45 => 0,
	50 => 0,
	60 => 0,
	70 => 0,
	80 => 0,
	90 => 0,
	100 =>0,
);

my %intraLenCountHash = %interLenCountHash;

my %interLenHash = (
  1 => 0,
  2 => 0,
  4 => 0,
  6 => 0,
  8 => 0,
  10 => 0,
  15 => 0,
  20 => 0,
  25 => 0,
  30 => 0,
  35 => 0,
  40 => 0,
  45 => 0,
  50 => 0,
  60 => 0,
  70 => 0,
  80 => 0,
  90 => 0,
  100 => 0,
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
	print "$nameHash{$k}\t$interLenCountHash{$k}\t$intraLenCountHash{$k}\t$total\t$interLenHash{$k}\t$intraLenHash{$k}\t$tlen\n";
}
#-- subs ---- 
sub sortLength{
	my($len, $type) = @_;
	if ($type eq "inter"){
		if ($len > 999 && $len <2000){
				$interLenCountHash{1}++;
				$interLenHash{1} += $len;
		}elsif ($len > 1999 && $len < 4000){
				$interLenCountHash{2}++;
				$interLenHash{2} += $len;
		}elsif ($len > 3999 && $len < 6000){
        $interLenCountHash{4}++;
				$interLenHash{4} += $len;	
		}elsif ($len > 5999 && $len < 8000){
        $interLenCountHash{6}++;
			$interLenHash{6} += $len;
		}elsif ($len > 7999 && $len < 10000){
        $interLenCountHash{8}++;
				$interLenHash{8} += $len;
		}elsif ($len > 9999 && $len < 15000){
        $interLenCountHash{10}++;
				$interLenHash{10} += $len;
		}elsif ($len > 14999 && $len < 20000){
        $interLenCountHash{15}++;
				$interLenHash{15} += $len;
		}elsif ($len > 19999 && $len < 25000){
        $interLenCountHash{20}++;
				$interLenHash{20} += $len;
		}elsif ($len > 24999 && $len < 30000){
        $interLenCountHash{25}++;
				$interLenHash{25} += $len;
		}elsif ($len > 29999 && $len < 35000){
        $interLenCountHash{30}++;
				$interLenHash{30} += $len;
		}elsif ($len > 34999 && $len < 40000){
        $interLenCountHash{35}++;
				$interLenHash{35} += $len;
		}elsif ($len > 39999 && $len < 45000){
        $interLenCountHash{40}++;
				$interLenHash{40} += $len;
		}elsif ($len > 44999 && $len < 50000){
        $interLenCountHash{45}++;
				$interLenHash{45} += $len;
		}elsif ($len > 49999 && $len < 60000){
        $interLenCountHash{50}++;
				$interLenHash{50} += $len;
		}elsif ($len > 59999 && $len < 70000){
        $interLenCountHash{60}++;
				$interLenHash{60} += $len;
		}elsif ($len > 69999 && $len < 80000){
        $interLenCountHash{70}++;
				$interLenHash{70} += $len;
		}elsif ($len > 79999 && $len < 90000){
        $interLenCountHash{80}++;
				$interLenHash{80} += $len;
		}elsif ($len > 89999 && $len < 100000){
        $interLenCountHash{90}++;
				$interLenHash{90} += $len;
		}elsif ($len > 99999 && $len < 4000000){
        $interLenCountHash{100}++;
				$interLenHash{100} += $len;
		}

	}elsif($type eq "intra"){
		if ($len > 999 && $len <2000){
				$intraLenCountHash{1}++;
				$intraLenHash{1} += $len;
		}elsif ($len > 1999 && $len < 4000){
				$intraLenCountHash{2}++;
				$intraLenHash{2} += $len;
		}elsif ($len > 3999 && $len < 6000){
        $intraLenCountHash{4}++;
				$intraLenHash{4} += $len;	
		}elsif ($len > 5999 && $len < 8000){
        $intraLenCountHash{6}++;
			$intraLenHash{6} += $len;
		}elsif ($len > 7999 && $len < 10000){
        $intraLenCountHash{8}++;
				$intraLenHash{8} += $len;
		}elsif ($len > 9999 && $len < 15000){
        $intraLenCountHash{10}++;
				$intraLenHash{10} += $len;
		}elsif ($len > 14999 && $len < 20000){
        $intraLenCountHash{15}++;
				$intraLenHash{15} += $len;
		}elsif ($len > 19999 && $len < 25000){
        $intraLenCountHash{20}++;
				$intraLenHash{20} += $len;
		}elsif ($len > 24999 && $len < 30000){
        $intraLenCountHash{25}++;
				$intraLenHash{25} += $len;
		}elsif ($len > 29999 && $len < 35000){
        $intraLenCountHash{30}++;
				$intraLenHash{30} += $len;
		}elsif ($len > 34999 && $len < 40000){
        $intraLenCountHash{35}++;
				$intraLenHash{35} += $len;
		}elsif ($len > 39999 && $len < 45000){
        $intraLenCountHash{40}++;
				$intraLenHash{40} += $len;
		}elsif ($len > 44999 && $len < 50000){
        $intraLenCountHash{45}++;
				$intraLenHash{45} += $len;
		}elsif ($len > 49999 && $len < 60000){
        $intraLenCountHash{50}++;
				$intraLenHash{50} += $len;
		}elsif ($len > 59999 && $len < 70000){
        $intraLenCountHash{60}++;
				$intraLenHash{60} += $len;
		}elsif ($len > 69999 && $len < 80000){
        $intraLenCountHash{70}++;
				$intraLenHash{70} += $len;
		}elsif ($len > 79999 && $len < 90000){
        $intraLenCountHash{80}++;
				$intraLenHash{80} += $len;
		}elsif ($len > 89999 && $len < 100000){
        $intraLenCountHash{90}++;
				$intraLenHash{90} += $len;
		}elsif ($len > 99999 && $len < 4000000){
        $intraLenCountHash{100}++;
				$intraLenHash{100} += $len;
		}
	}
}
