package Blasted;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use Exporter;
$VERSION=0.50;
@ISA=qw(Exporter);

@EXPORT =qw();			#use Blast;
@EXPORT_OK = qw(&parse_query);	#use Blast (...);
%EXPORT_TAGS = ();   #use Blast qw(:TAG1);

####################################
##########CODE #####################

sub parse_query { # FILEHANDLE 
 	my %args = (FILEHANDLE => '', ALIGN => 0, @_ ); #FH is glob filehadnele \*IN
 	die "Expecting a glob to the filehandle e.g. FH => \\*IN\n" if $args{'FILEHANDLE'} !~ /GLOB/;
	my $FH = $args{'FILEHANDLE'};  #too long type all the time
	my %d=();
	my $line = '';
	$line = <$FH> until $line=~/^(BLAST[A-Z.0-9]+)/ || eof $FH;
	if (eof $FH ) { 
		$d{'error'} = 'No Blast Title Found';
		return \%d;
	}
	$d{'program'} = $1;

	$line =<$FH> until $line=~ /Query= ([0-9A-Za-z_.]+)/;
	$d{'qname'}=$1;
	
	my $tmp='';
	until ( $line =~ /letters\)/ ) {
		$tmp.=$line;
		$line =<$FH>;
	}
	$line =~ /([0-9,]+)/;
	$d{'qlen'} =$1;
	$d{'qlen'} =~ s/\,//m;
	$tmp=~s/Query= +$d{'qname'} *//;
	$tmp=~s/\n/ /mg;
	$tmp=~s/ +/ /mg;
	$tmp=~s/^\|//;
	$tmp=~s/ +$|^ +//m;
	$d{'qdefn'}=$tmp;
	$line=<$FH> until $line =~/^Database: +(\S+)/;
	$d{'db'}{'name'} = $1;
	$line=<$FH>;
	($d{'db'}{'seq#'}, $d{'db'}{'len'}) = $line=~ /([0-9,]+) +seq.* ([0-9,]+)/;
	$d{'db'}{'seq#'} =~ s/\,//m;
	$d{'db'}{'len'} =~ s/\,//m;
	my $s=0;
	my @sbjct;
	until ( $line=~/^\s+Database/ ) {
		$line=<$FH> until $line=~/^(>)/ || $line=~/^\s+Database/;
 		next if $1 ne '>';
 		#we hit a >subject
 		($sbjct[$s]{'name'} )= $line=~/>([0-9A-Za-z_.]+)/;
 		my $tmp='';
		until ( $line =~ /Length = / ) {
			$tmp.= $line;
			$line =<$FH>;
		}
		$tmp=~s/\n/ /mg;
		$tmp=~s/ +/ /mg;
		$tmp=~s/>$sbjct[$s]{'name'} *//;
		$tmp=~s/^\|//;
		$tmp=~s/ +$|^ +//m;
		$sbjct[$s]{'defn'}=$tmp;
		
		
		$line=~ /= ([0-9,]+)/;
		$sbjct[$s]{'len'}=$1;
		$sbjct[$s]{'len'}=~s/\,//m;
		my @hsp;
		my $h=0;
		until ( $line=~ /^>/ || $line=~ /\s+Database/ ) {
			$line=<$FH> until $line=~ /^\s+(Score) =/ || $line=~ /^>/ || $line=~ /\s+Database/; 
			next if $1 ne 'Score';
			( $hsp[$h]{'bits'}, $hsp[$h]{'score'} ) = $line=~ /Score = (\S+) bits \((\d+)/;
			( $hsp[$h]{'expect'}) = $line =~ /Expect = (\S+)/;
			my $tmp = <$FH>;
			$line = <$FH>;
			$tmp.= $line;   #wordwrap occurs in big files with ncbi blast
			$tmp=~s/\n/ /mg;	
			($hsp[$h]{'bpident'},$hsp[$h]{'bpalign'}, $hsp[$h]{'%ident'}) 
					= $tmp=~/Identities = (\d+)\/(\d+) +\(([0-9.]+)/;
			($hsp[$h]{'bppos'},$hsp[$h]{'%pos'} ) = $tmp=~/Positivies = (\d+)\/\d+ +\(([0-9.]+)/;
			($hsp[$h]{'bpgap'}, $hsp[$h]{'%gap'} ) = $tmp=~/Gaps = (\d+)\/\d+ +\(([0-9.]+)/;
			##parse pairwise alignment ##
			my (@QA, @SA, @MA);
			my ($qb,$qe,$sb,$se) = (0,0,0,0);
			$line=<$FH> until $line =~/Query:/;  #jump to fist alignment line
			until ($line !~ /^Query:/) {
				$line =~ /(Query:\s+)(\d+)\s+(\S+)\s+(\d+)/;
				$qb=$2 unless $qb;   push @QA, $3; $qe=$4;
				$line=<$FH>;
				push @MA, substr($line, length($1)-1,length($3) );
				$line=<$FH>;
				$line=~/Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/;
				$sb=$1 unless $sb;	push @SA, $2; $se=$3;
				$line=<$FH>;
				$line=<$FH> until $line=~/\S/;
			}
			$hsp[$h]{'qb'}=$qb;  $hsp[$h]{'qe'}=$qe;
			$hsp[$h]{'sb'}=$sb;	$hsp[$h]{'se'}=$se;
			if ($args{'ALIGN'} ) {
				$hsp[$h]{'alignq'}=join ("", @QA);
				$hsp[$h]{'alignm'}=join ("",@MA);
				$hsp[$h]{'aligns'}=join ("",@SA);
			}
				
			$h++;
		}
		@{$sbjct[$s]{'hsps'} }=@hsp;
		$s++;
	}
	@{ $d{'sbjct'} }=@sbjct;
 	$line =<$FH> until $line=~/^BLAST/ || eof $FH;
  	return \%d;
}
	






######last line #####
1;
