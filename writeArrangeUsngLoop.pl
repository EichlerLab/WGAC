

#!/usr/bin/perl
# this is used to write arrange fil


if (($ARGV[0] eq "-h") || (!$ARGV[0])){
	print " this is used to generate a arranged file ";
	print " using showseq file. \n";
	print " Usage: writeArrangeUsngLoop.pl showseq.out\n";
	exit;
} 

my $showseq=$ARGV[0];
open(SF,$showseq) || die "Can not open $showseq\n";

while(<SF>){
	if(/seq/ || /SEQ/){
		next;
	}
	chomp; 
	my @line=split(/\t/, $_);
	my $chr=$line[0];
	print "NEWLINE\n$chr\t1\n";
}
close SF;
