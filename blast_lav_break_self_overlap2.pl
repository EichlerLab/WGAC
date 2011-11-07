#!/usr/bin/perl
#THIS PROGRAM PARSES BLASTOUTPUT IN A HIT BY HIT BASIS
#00-04-15 added "}" printing if no alignments #
#00-04-15 added eof checking in case of no alignments#
use strict 'vars';

use Getopt::Long;
use Data::Dumper;
use lib '/home/xws/bin/JABPerlMod';
use Blast qw(&parse_query);

use vars qw($true $false);
$true=1; $false=0;


use vars qw(%opt);
use vars qw(%prev_queries @files $path @defaults);


if ($ARGV[0] eq '') {
	print "blast_lav_break_self_overlap.pl
  --in [file/directory] containing blast records to parse
  --out [directory] files will be *.brk)
  --options
  ";
	die "";
}
####GET OPTIONS
&GetOptions(\%opt, "in=s", "out:s","options:s");
$opt{'in'} || die "Please designate input file with -in\n";
$opt{'options'} || ( $opt{'options'} = "" ) 	;
$opt{'options'} = ', '.$opt{'options'} if $opt{'options' ne ''};
$opt{'out'} || ($opt{'out'} = '.');

 
@defaults= ( );
#		'MAX_%GAP'=> 80, MIN_BPALIGN => 10, MIN_FRACBPMATCH=>0.80,
#		,SKIP_OVERLAP=> $false, SKIP_SELF => $false, SKIP_IDENT_HSP=> $false) };
push @defaults, split " *[=>,]+ *",$opt{'options'}  if $opt{'options'};

%prev_queries=();
if (opendir (DIR, $opt{'in'}) ) {
	@files = grep { /[a-zA-Z0-9]/ } readdir DIR;
	close DIR;
	$path = $opt{'in'};
} elsif (open (IN, $opt{'in'})) {
	($path ) = $opt{'in'} =~ /(^.*)\//;
	$opt{'in'} =~ s/^.*\///;
	@files=($opt{'in'});
} else {
	die "--in  ($opt{'in'})  not a file and not a directory\n";
}

foreach my $f (@files) {
	print "$path/$f\n";
	open (IN, "$path/$f") || die "Can't open file ($f)! \n";
	open (OUT, ">$opt{'out'}/$f.brk") || die "Can't open file ($f.brk)!\n";
	&lav_remove_self_overlap( FILEHANDLE => \*IN , OUTHANDLE => \*OUT,@defaults);	
	close IN;
	close OUT;
}

sub lav_remove_self_overlap{ # FILEHANDLE 
 	my %args = (FILEHANDLE => '', #FH is glob filehandle \*IN
 					OUTHANDLE => '', #OH is glob filehandel \*OUT
 		@_ ); 
	die "Expecting a glob to the filehandle e.g. FH => \\*IN\n" if $args{'FILEHANDLE'} !~ /GLOB/;
	die "Expecting a glob to the outfile handle e.g.OUTHANDLE => \\*IN\n" if $args{'OUTHANDLE'} !~ /GLOB/;

	my $FH = $args{'FILEHANDLE'};  #too long type all the time
	my $OH = $args{'OUTHANDLE'};

	my $line = '';	
	$line = <$FH> until $line=~ /^\#\:lav/ || eof $FH;
   die "End of file reached without finding #:lav header\n" if (eof $FH);
	#START PARSING SINGLE QUERY VS SINGLE SUBJECT WITH POSSIBLE REVERSE SUBJECT#
	my $search_string='';
	until ( $line=~/^s {/) {
		print $OH $line;
		$search_string.=$line;
		$line =<$FH>;
		
	}
	my ($m,$i,$v,$o,$e)=$search_string=~/M = (\d+), I = (-\d+), V = (-\d+), O = (\d+), E = (\d+)/;
	print "SCORES: $m, $i, $v, $o, $e\n";
	
	print $OH $line;
	$line =<$FH>;
	my $line2 = <$FH>;
	print $OH $line, $line2;
	if ($line eq $line2) {
	#####same piece same orientation, so break it up!	
		print "We have found self!\n";
		$line = <$FH>;
		until ($line=~ /^a \{/ || eof ($FH) ) {
			last if eof $FH;
			print $OH $line;
			$line =<$FH>;
		}
		print $OH "}" if eof ($FH);  #because if eof line will have nothing
		
		
		while ($ line =~ /^a \{/) {
			#print "START:$line";
			my @l =();
			push @l, $line;
			$line =<$FH>;
			until ( $line =~ /^\w \{/  || $line =~ /\#\:lav/ || eof ($FH) ) {
				push @l, $line;
				$line =<$FH>;
			}
			#lets calculate#
			my ($s) = $l[1] =~ /s (\d+)/;
			my ($b1,$b2)= $l[2]=~ /b (\d+) (\d+)/;
			my ($e1,$e2)= $l[3]=~ /e (\d+) (\d+)/;
			
			#remove those extra lines#

				
			if ($b2 <=$e1) {
				shift @l; shift @l; shift @l; shift @l; pop @l;
				print "BAD\n";
				my @fpieces;
				print "***NEW ALIGNMENT $s $b1($b2)  $e1($e2)  ";
				LINES: while (@l > 0 ) {
					my @piece=();
					my $starting= $true;
					my $begin2=-99;
					foreach my $l (@l) {
						next if $l !~ /l (\d+) (\d+) (\d+) (\d+) (\d+)/;
						my ($qb,$sb,$qe,$se,$perc)=($1,$2,$3,$4,$5);
						if ($starting==$true) {
							$begin2 = $sb;
							$starting=$false;
							print scalar(@l), "PIECE IS \n";
						}
						#check if HSP overlapping#
						last if ($qe >= $sb);
						last if ($qe >= $begin2);
						push @piece, $l;
						#print "$l";
					}
					push @fpieces, \@piece if @piece > 0;
					shift @l;
				}

				##################loop to score and print####
				foreach my $f (@fpieces) {
					my @piece = @{$f};
					my ($bpident,$bpgap,$bpalign, $gaps)=(0,0,0,0);
					my ($lqe,$lse)=(0,0);
					my ($q0,$s0, $qn, $sn)=(0,0);
					foreach my $p (@piece) {
						next if $p !~ /l (\d+) (\d+) (\d+) (\d+) (\d+)/;
						my ($qb,$sb,$qe,$se,$perc)=($1,$2,$3,$4,$5);
						$q0=$qb if $q0==0;
						$s0=$sb if $s0==0;
						($qn,$sn)=($qe,$se);
						my $leng= $qe-$qb+1;
						##this is the best guess since he is doing alot of rounding##
						#print "$line\n";
						if ($perc==0 || $leng<1) {
							print "BAD ALIGNMENT$p ($leng)\n";
							next;
						}
						my $ident=int ($leng*($perc+0.499999)/100 ); #reverse engineering#
						my $fix=int($ident/$leng *100+ 0.5);  #forward check#
						#son of a bitch#
						#print "$leng ($perc) => $ident  fix ($fix)\n";
						$bpident+=$ident;
						$bpalign+=$leng;
						if ($lqe>0) {
							###calculate gap statistics
							my $qgap=$qb-$lqe-1;
							my $sgap=$sb-$lse-1;
							if ($qgap>0 and $sgap>0) {
								print "PRINT DOUBLE GAP ERROR\n" ;
							}
							$bpgap+=$qgap + $sgap;
							$gaps++;
						}
						($lqe,$lse)=($qe,$se);
					}
					####print scores
					my $newscore=$bpident *$m+($bpalign-$bpident)*$v-$bpgap*$e-$o*$gaps;
					#print "===>$bpident, $bpalign, $bpgap ($gaps) NEWSCORE=$newscore\n";
					print $OH "a {      #broken $s\n";
					print $OH "  s $newscore\n";
					print $OH "  b $q0 $s0\n";
					print $OH "  e $qn $sn\n";
					foreach (@piece) { print $OH $_;}
					print $OH "}\n";
		
				}
		
			} else {
				#######GOOD NONOVERLAPPING FORWARD ALIGNMENT#####
				foreach my $l (@l) { print $OH $l;}
			}

		}
	}
	
	#REVERSE ORIENTATION AOK#
	while (<$FH>) {print $OH $_;}
	return;

} 

