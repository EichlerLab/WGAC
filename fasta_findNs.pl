#!/usr/bin/perl
#LOAD MODULES
#00-04-15 This program begs to be cleaned up.
use Getopt::Std;
use strict 'vars';




use vars qw($true $false);
$true=1; $false=0;
use vars qw($opt_i $opt_p $opt_o $opt_s);
use vars qw(@fastas);

if ($ARGV[0] eq '') {

	print "fasta_findNs.pl ******************************************\n";
	print "  This program extracts a subsequence from a fasta record.\n";
	print "\t-i [path] fasta file or directory of fasta records\n";
	print "\t          or : delimited files\n";
	print "\t	   or the first column of a list\n";
	print "\t-p [path] directory to find fasta files in (optional)\n";
	print "\t-o [path] output table (default Npositions.tbl)\n";
	print "\t-s [integer] number of rows to initially skip in list (default zero)\n";
		die "***********************************************************\n";
}
getopts('i:p:s:o:');
$opt_i || die "Please enter path for directory,fasta(s), or list";
$opt_p ||= '';
$opt_s ||=  0;
$opt_o ||= 'Npositions.tbl';
############################
if (opendir (DIR, $opt_i)) {
	print "OPENING DIRECTORY $opt_i\n";
	@fastas = grep {/[a-zA-Z0-9]/} readdir(DIR);
	$opt_p=$opt_i;
} elsif ($opt_i =~ /\:/) {
	print "SPLITING MULTIPLE FILE NAMES\n";
	@fastas=split ":", $opt_i;
} else {
	open (TEST,$opt_i) || die "Can not open supposed table, fasta or directory $opt_i\n";
	my $header=<TEST>;
	print "SINGLE FASTA RECORD";
	if ($header=~/^>/) {
		#single fasta record
		$fastas[0]=$opt_i;
	} else {
		print "ASSUMING THAT RECORD IS TABLE\n";
		my $line=$header;
		#assuming this is a table
		for (my $i=0; $i<$opt_s; $i++) { $line=<TEST>;}
		$line=~ s/\n/\t/;
		my @col=split "\t",$line;
		push @fastas, $col[0] if $col[0] != 0;
	}
}
$opt_p.='/' if $opt_p;
print "OPEN $opt_o\n";
open(OUT,">$opt_o") || die "Can't open output file $opt_o\n";
print OUT "NAME\tbegin\tend\tcolor\twidth\tdefn\n";

foreach my $f (@fastas) {
	open(IN, $opt_p.$f ) || die "Can't open $opt_p" . "$f\n";
	print "SEARCHING FOR Ns in $f....\n";
	my @ns=&fasta_ns_positions(-filehandle=>\*IN,-min_size=>10);
	$f=~s/^.*\///;
	foreach my $n (@ns) {
	
		print OUT "$f\t$$n{'begin'}\t$$n{'end'}\tblack\t7\tpolyNtract\n";
	}
	
}
	



sub fasta_ns_positions {
	my %args=(-filehandle=>\*STDIN, -min_size=>1 ,@_);
	my @ns=();
	my $position=0;
	my $i=0;
	my $fh=$args{'-filehandle'};
	my $min_size = $args{'-min_size'};
	my $header= <$fh>;
	my $new_begin=$true;
	while ( <$fh> ) {
		s/\r\n/\n/;
		chomp; chomp;
		s/\s+//;
		my $llen=length;
		$position += $llen;
		my $line = $_;
		#print "$line  P:$position =>";
		my $ns_in_line=$false;
		while ( $line=~ /([nN]+)/g ) {
			$ns_in_line=$true;
			#print "FOUND",length($1),"  ",pos $line,"";
			if ( ( pos($line)-length($1) ) > 0  ||  ((length($1)- pos $line )==0 && $new_begin==$true) ) {
				#print "BEGIN\n";
				$ns[$i]{'begin'}=$position-$llen+(pos $line)-(length($1))+1;
				$new_begin=$false;
			}
			if (pos $line < length($line)) {
				$ns[$i]{'end'}=$position-$llen+pos $line;
				#print "\n$i: $ns[$i]{'begin'} $ns[$i]{'end'} ",($ns[$i]{'end'}-$ns[$i]{'begin'}+1);
				$i++ if ($ns[$i]{'end'}-$ns[$i]{'begin'}+1) >= $min_size;
				$new_begin=$true;
				#print " $i";
				#my $pause=<STDIN>;
			}
			
		}
		if ($new_begin==$false) { $ns[$i]{'end'}=$position;}
		if ($new_begin==$false && $ns_in_line==$false) { 
			$ns[$i]{'end'}=$position-$llen;
			#print "\n$i: $ns[$i]{'begin'} $ns[$i]{'end'} ",($ns[$i]{'end'}-$ns[$i]{'begin'}+1);
			$i++ if ($ns[$i]{'end'}-$ns[$i]{'begin'}+1) >= $min_size;
			$new_begin=$true;
			#print " $i";
			#my $pause=<STDIN>;
		}
		#print "\n";
	}
	return @ns;
	
}

