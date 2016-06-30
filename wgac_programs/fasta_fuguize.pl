#!/usr/bin/perl
#01-10-17 fixed fasta_save to include the last base on a word wrap
#00-06-01 renmae fasta_fuguize.pl
#00-05-21 removed references to unused columns in repeatmasker table#

#!/usr/local/bin/perl
#LOAD MODULES

use Getopt::Std;
use strict 'vars';

use vars qw($true $false);
$true=1; $false=0;

use vars qw($opt_f $opt_r $opt_o $opt_g);
use vars qw($title $position $subseq $header);

if ($ARGV[0] eq '') {

	print "fasta_fuguize.pl ******************************************\n";
	print "This program encodes/fuguizes a fasta file given repeatmasker outfile.\n";
	print "-f [path] fasta sequence file is required\n";
	print "-r [path] repeat masker *.out file\n";
	print "-o [path] output for encoded file\n";
	print "-g [switch] requires only -f and adds .out .fugu\n";
	print "*This program does not change the case of sequence\n";
	print "***********************************************************\n";
	exit;
}

getopts('f:r:o:g');
if ($opt_g) {
	$opt_r= $opt_f. '.out';
	$opt_o= $opt_f. '.fugu';
}
my ($fasta,$fasta_header) = &load_fasta($opt_f);
my ($fastalength)=length ($fasta);
chomp $fasta_header;
chomp $fasta_header;
$fasta_header .= " UELEN:$fastalength";
print "FASTA length ", $fastalength,"===>";
my ($efasta) = &fasta_encode($fasta,$opt_r,$opt_o);
print "FUGU  length ", length($efasta),"\n";
&save_fasta("$opt_o",$efasta, $fasta_header);


####decode test ####
#my ($b,$e) = &rm_out_decode(5,30437, $opt_r);
#print "$b to $e \n";

#####SUBROUTINES ###############
sub fasta_encode {
	my ($fasta, $rm_out, $fasta_out) = @_;
	my @repeats= &load_rm_out($rm_out);
	@repeats =reverse @repeats;
	
	foreach my $r (@repeats) {
		substr( $fasta, $$r{'b'}-1, $$r{'e'}- $$r{'b'}+1 ) ='';
		#print "$$r{'b'} $$r{'e'} $$r{'rep'} $$r{'fam'}\n";
	}
	return $fasta	
}

		
sub rm_out_decode {
	my ($begin, $end, $rm_out) = @_;
	my @repeats = &load_rm_out($rm_out);
	###need more code to do this
	foreach my $r (@repeats) {
		my $l= $$r{'e'}-$$r{'b'} +1;
		$begin += $l if $begin >= $$r{'b'} ;
		 if ($end >= $$r{'b'}) {
				$end += $l;
		} else {
			last;
		}
	}
				
	return ($begin, $end);
}

sub load_fasta {
	open (FASTAIN, "$_[0]") || die "Can't open $_[0]\n";
	my $fasta = '';
	my $header = <FASTAIN>;
	while (<FASTAIN>) {
		s/\r\n/\n/;
		chomp;
		$fasta .= $_;
	}
	return ($fasta, $header);
	
}
	
sub load_rm_out {
		open (RMOUT, "$_[0]")  || die "Can't open $_[0]\n";
		my @repeats=();
		my $line =<RMOUT>; $line=<RMOUT>; $line=<RMOUT>;
		my $i;
		while (<RMOUT>) {
			###need to read the file backward
			chomp;
			s/^ +//;
			my @col=split " +";
			if ($i>0 && ($col[5] <= $repeats[$i-1]{'e'} ) ) {
				#print "OVERLAP\n";
				$i--;
			} else {
				$repeats[$i]{'b'} = $col[5];
			}
		
			$repeats[$i]{'e'} = $col[6];
			#not used in this implementation#
			#$repeats[$i]{'rep'} .= "_$col[9]";
			#$repeats[$i]{'fam'} .= "_$col[10]";
			#print $repeats[$i]{'b'}, "\n";
			$i++;
		}
		close RMOUT;
	   return @repeats;
}

sub save_fasta {
	my ($filename, $whole_seq, $header) = @_;
	open (FASTAOUT, ">$_[0]") || die "Can't create $_[0]\n";
	print FASTAOUT "$header\n";

	my $width=60;
	my $m=0;
	for ($m=0; $m+$width<length ($whole_seq); $m+=$width) {
		print  FASTAOUT substr($whole_seq, $m, $width),"\n";
	}
	if ($m <=length($whole_seq)-1) {
		print FASTAOUT  substr($whole_seq, $m),"\n";
	}
	close FASTAOUT;
}





