

#!/usr/bin/perl
use strict;
my @a;
my $intraCount = 0;
my $interCount = 0;
my ($intraSim, $interSim, $plen);


## first to set up variables
my %interSimCountHash = (
	0.90 => 0,
	0.905 => 0,
	0.91 => 0,
	0.915 => 0,
	0.92 => 0,
	0.925 => 0,
	0.93 => 0,
	0.935 => 0,
	0.94 => 0,
	0.945 => 0,
	0.95 => 0,
	0.955 => 0,
	0.96 => 0,
	0.965 => 0,
	0.97 => 0,
	0.975 => 0,
	0.98 => 0,
	0.985 => 0,
	0.99 => 0,
	0.995 => 0,
	1 => 0
);

my %intraSimCountHash = %interSimCountHash;

my %interSimLengthHash = (
	0.90 => 0,
	0.905 => 0,
	0.91 => 0,
	0.915 => 0,
	0.92 => 0,
	0.925 => 0,
	0.93 => 0,
	0.935 => 0,
	0.94 => 0,
	0.945 => 0,
	0.95 => 0,
	0.955 => 0,
	0.96 => 0,
	0.965 => 0,
	0.97 => 0,
	0.975 => 0,
	0.98 => 0,
	0.985 => 0,
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
print "size\tinter\tintra\tTotal\tinterLength\tintraLength\tTotalLength\n";
foreach $k(@keyArray){
	$total = $interSimCountHash{$k} + $intraSimCountHash{$k};
	$tlen = $interSimLengthHash{$k} + $intraSimLengthHash{$k};
    my $fPerCent = $k * 100;
	print "$fPerCent\t$interSimCountHash{$k}\t$intraSimCountHash{$k}\t$total\t$interSimLengthHash{$k}\t$intraSimLengthHash{$k}\t$tlen\n";
}
#-- subs ---- 
sub sortSim{
	my($sim, $type, $len) = @_;
	if ($type eq "inter"){
		if ($sim >= 0.9 && $sim < 0.905){
				$interSimCountHash{0.90}++;
				$interSimLengthHash{0.90} += $len;
		}elsif ($sim >= 0.905 && $sim < 0.91){
				$interSimCountHash{0.905}++;
				$interSimLengthHash{0.905} += $len;
		
		}elsif ($sim >= 0.91 && $sim < 0.915){
				$interSimCountHash{0.91}++;
				$interSimLengthHash{0.91} += $len;
		}elsif ($sim >= 0.915 && $sim < 0.92){
				$interSimCountHash{0.915}++;
				$interSimLengthHash{0.915} += $len;
		
		}elsif ($sim >= 0.92 && $sim < 0.925){
        $interSimCountHash{0.92}++;
				$interSimLengthHash{0.92} += $len;
		}elsif ($sim >= 0.925 && $sim < 0.93){
        $interSimCountHash{0.925}++;
				$interSimLengthHash{0.925} += $len;
		
		}elsif ($sim >= 0.93 && $sim < 0.935){
        $interSimCountHash{0.93}++;
				$interSimLengthHash{0.93} += $len;
		}elsif ($sim >= 0.935 && $sim < 0.94){
        $interSimCountHash{0.935}++;
				$interSimLengthHash{0.935} += $len;
		
		}elsif ($sim >= 0.94 && $sim < 0.945){
        $interSimCountHash{0.94}++;
				$interSimLengthHash{0.94} += $len;
		}elsif ($sim >= 0.945 && $sim < 0.95){
        $interSimCountHash{0.945}++;
				$interSimLengthHash{0.945} += $len;
		
		}elsif ($sim >= 0.95 && $sim < 0.955){
        $interSimCountHash{0.95}++;
				$interSimLengthHash{0.95} += $len;
		}elsif ($sim >= 0.955 && $sim < 0.96){
        $interSimCountHash{0.955}++;
				$interSimLengthHash{0.955} += $len;
		
		}elsif ($sim >= 0.96 && $sim < 0.965){
        $interSimCountHash{0.96}++;
				$interSimLengthHash{0.96} += $len;
		}elsif ($sim >= 0.965 && $sim < 0.97){
        $interSimCountHash{0.965}++;
				$interSimLengthHash{0.965} += $len;
		
		}elsif ($sim >= 0.97 && $sim < 0.975){
        $interSimCountHash{0.97}++;
				$interSimLengthHash{0.97} += $len;
		}elsif ($sim >= 0.975 && $sim < 0.98){
        $interSimCountHash{0.975}++;
				$interSimLengthHash{0.975} += $len;

		}elsif ($sim >= 0.98 && $sim < 0.985){
        $interSimCountHash{0.98}++;
				$interSimLengthHash{0.98} += $len;
		}elsif ($sim >= 0.985 && $sim < 0.99){
        $interSimCountHash{0.985}++;
				$interSimLengthHash{0.985} += $len;

		}elsif ($sim >= 0.99 && $sim < 0.995){
        $interSimCountHash{0.99}++;
				$interSimLengthHash{0.99} += $len;
		}elsif ($sim >= 0.995 && $sim <= 1){
				if ($sim == 1){
				$interSimCountHash{1}++;
				$interSimLengthHash{1} += $len;
		}else{
        $interSimCountHash{0.995}++;
				$interSimLengthHash{0.995} += $len;
				}
		}
	}elsif($type eq "intra"){
		if ($sim >= 0.9 && $sim < 0.905){
				$intraSimCountHash{0.90}++;
				$intraSimLengthHash{0.90} += $len;
		}elsif ($sim >= 0.905 && $sim < 0.91){
				$intraSimCountHash{0.905}++;
				$intraSimLengthHash{0.905} += $len;
		
		}elsif ($sim >= 0.91 && $sim < 0.915){
				$intraSimCountHash{0.91}++;
				$intraSimLengthHash{0.91} += $len;
		}elsif ($sim >= 0.915 && $sim < 0.92){
				$intraSimCountHash{0.915}++;
				$intraSimLengthHash{0.915} += $len;
		
		}elsif ($sim >= 0.92 && $sim < 0.925){
        $intraSimCountHash{0.92}++;
				$intraSimLengthHash{0.92} += $len;
		}elsif ($sim >= 0.925 && $sim < 0.93){
        $intraSimCountHash{0.925}++;
				$intraSimLengthHash{0.925} += $len;
		
		}elsif ($sim >= 0.93 && $sim < 0.935){
        $intraSimCountHash{0.93}++;
				$intraSimLengthHash{0.93} += $len;
		}elsif ($sim >= 0.935 && $sim < 0.94){
        $intraSimCountHash{0.935}++;
				$intraSimLengthHash{0.935} += $len;
		
		}elsif ($sim >= 0.94 && $sim < 0.945){
        $intraSimCountHash{0.94}++;
				$intraSimLengthHash{0.94} += $len;
		}elsif ($sim >= 0.945 && $sim < 0.95){
        $intraSimCountHash{0.945}++;
				$intraSimLengthHash{0.945} += $len;
		
		}elsif ($sim >= 0.95 && $sim < 0.955){
        $intraSimCountHash{0.95}++;
				$intraSimLengthHash{0.95} += $len;
		}elsif ($sim >= 0.955 && $sim < 0.96){
        $intraSimCountHash{0.955}++;
				$intraSimLengthHash{0.955} += $len;
		
		}elsif ($sim >= 0.96 && $sim < 0.965){
        $intraSimCountHash{0.96}++;
				$intraSimLengthHash{0.96} += $len;
		}elsif ($sim >= 0.965 && $sim < 0.97){
        $intraSimCountHash{0.965}++;
				$intraSimLengthHash{0.965} += $len;
		
		}elsif ($sim >= 0.97 && $sim < 0.975){
        $intraSimCountHash{0.97}++;
				$intraSimLengthHash{0.97} += $len;
		}elsif ($sim >= 0.975 && $sim < 0.98){
        $intraSimCountHash{0.975}++;
				$intraSimLengthHash{0.975} += $len;

		}elsif ($sim >= 0.98 && $sim < 0.985){
        $intraSimCountHash{0.98}++;
				$intraSimLengthHash{0.98} += $len;
		}elsif ($sim >= 0.985 && $sim < 0.99){
        $intraSimCountHash{0.985}++;
				$intraSimLengthHash{0.985} += $len;

		}elsif ($sim >= 0.99 && $sim < 0.995){
        $intraSimCountHash{0.99}++;
				$intraSimLengthHash{0.99} += $len;
		}elsif ($sim >= 0.995 && $sim <= 1){
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
