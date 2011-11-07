

#!/usr/bin/perl

## this is used to split a pair 
## lable the inter, intra
## and to switch the start and end if they are on the reverse strain
## so that they are going to be sorted correctly

use strict;
my @a;
my $id = "id";  # pair id
my $count = 0;
my $qid;
my $sid;
my $ori1 = "+";
my $ori2 = "+";

while(<>){
	next if (/^\n/ || /NAME/);
  if(/QNAME/){
			print "NAME\tBegin\tEnd\tChromosomeLen\tLength\tSimilarity\tPairid\tPairtype\n";
			next;
	}    
	@a = split (/\t/);
			$count++;
			if ($a[1] > $a[2]){	# change start and end
					($a[1], $a[2])=	switchOrder($a[1], $a[2]);
					$ori1 = "-";
			}else{
					$ori1 = "+";
			}

			if ($a[5] > $a[6]){	# change start and end
					($a[5], $a[6])= switchOrder($a[5], $a[6]);
					$ori2 = "-";
     			}else{
					$ori2 = "+";
			}
# print "$a[0] and $a[4]-------------\n";
			
			$id = "p$count";
			if ($a[0] eq $a[4]){	
					print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[22]\t$a[27]\t$id\tintra\t$ori1\n";
					print "$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[22]\t$a[27]\t$id\tintra\t$ori2\n";
			}else{
						$qid = $a[0];
						$sid = $a[4];
						$qid =~ s/_random//;
						$sid =~ s/_random//;
						if ($qid eq $sid){
								print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[22]\t$a[27]\t$id\tintra\t$ori1\n";
			          				print "$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[22]\t$a[27]\t$id\tintra\t$ori2\n";
						}else{
								print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[22]\t$a[27]\t$id\tinter\t$ori1\n";
      								print "$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[22]\t$a[27]\t$id\tinter\t$ori2\n";
						}
			}
}


sub switchOrder{
	my ($b, $e) = @_;
	my $tmp = $b;
#	print " ================= $b and $e \n";
	$b = $e;
	$e = $tmp;
	return($b, $e);
}

			
