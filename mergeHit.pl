#!/usr/bin/perl
## This is used to merge hits on file one line per hit
## it has to work on the presort file

use strict;
my @a;
my @end;
my @begin;
my $chrlen = 0;
my $pid = "chrid";

while(<>){

	next if (/NAME/ || /\=/);
	chomp;
	@a = split (/\s+/);
		if ($a[0] eq $pid){
						# continue to collect the begin and end
					push(@begin,$a[1]);
					push(@end, $a[2]);
		}else{
					if($a[3]){
						processNprint();
						$pid = $a[0];
						undef @begin;
						undef @end;
						push(@begin,$a[1]);
          					push(@end, $a[2]);
						$chrlen = $a[3];
					}
		}
}
processNprint();
sub processNprint{
	## the @begin and @end should have the saem number of coordinates
	my (@fbegin, @fend);
	my $begin = $begin[0];
	my $end = $end[0];		## to get the first pair
	my $arrayLen = @begin;
$arrayLen++;
	my $tmpN = 0;
	my $i= 1;
	my $mergeLen = 0;

	FORLOOP:for($i = 1; $i <=  $arrayLen+1; $i++){		# start with 1 because the first has taken
				if (!$begin[$i]){
						$mergeLen = abs ($begin - $end);
						 print "$pid\t$begin\t$end\t$mergeLen\t$chrlen\n";
					#	print "$begin[$i-1] and $arrayLen+1\n";
						last FORLOOP;
				}
				if ($begin[$i] > $end[$i]){			# change the begin and end
						$tmpN = $begin[$i];
						$begin[$i] = $end[$i];
						$end[$i] = $tmpN;
				}

				if($begin[$i] <= $end){			# overlap , merge
							if($end[$i] > $end){
										$end = $end[$i];
							}else{
										# still keep the the $end
							}
				}else{
						## now they are not overlap
						$mergeLen = $end - $begin;
						print "$pid\t$begin\t$end\t$mergeLen\t$chrlen\n";
						$begin = $begin[$i];
						$end = $end[$i];
				}
	}
}



