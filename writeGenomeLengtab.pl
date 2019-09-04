#!/usr/bin/perl
my @a;
my $hashRef;
my $startArrayRef;
my $endArrayRef;
my $pa2 = 0;
my $pa3 = 0;

print "seqname\tlength\n";
while(<>){
next if (/TOTAL/);
	s/\.fa//;
	@a = split (/\s+/);
	print "$a[1]\t$a[3]\n";
}



