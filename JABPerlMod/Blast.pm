package Blast;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

#00-05-27  added MIN and MAX HSPS #######

use Exporter;
$VERSION=0.50;
@ISA=qw(Exporter);

@EXPORT =qw();			#use Blast;
@EXPORT_OK = qw(&parse_query);	#use Blast (...);
%EXPORT_TAGS = ();   #use Blast qw(:TAG1);

####################################
##########CODE #####################
use vars qw($true $false);
($true,$false)=(1,0);

sub parse_query { # FILEHANDLE 
 	my @starts = (FILEHANDLE => '', #FH is glob filehandle \*IN
 				ALIGN => $false, SKIP_OVERLAP => $false, 
 				SKIP_ALLOVERLAP_ORIENT => $false,
 				SKIP_IDENT_HSP=> $false, SKIP_SUBJECTS=> '', SKIP_SELF=> $false,
 				
 				MIN_SIZEALIGN=> 0,  MAX_SIZEALIGN => 9999999999,
 				MIN_BPIDENT => 0,   MAX_BPIDENT =>  9999999999,
 				'MIN_%IDENT' => 0,  'MAX_%IDENT' => 101,
 				MIN_BPPOS => 0,  MAX_BPPOS => 9999999999,
 				'MIN_%POS' => 0, 'MAX_%POS' =>101,
 				MIN_BPGAP =>0, MAX_BPGAP => 999999999999,
 				'MIN_%GAP' => 0, 'MAX_%GAP' =>101,
 				
 				MIN_BPALIGN=> 0,  MAX_BPALIGN => 9999999999,    # matches plus mismatches 
 				MIN_FRACBPMATCH => 0, MAX_FRACBPMATCH =>1.01,   # matches / (matches + mismatches)	
 				
 				MIN_SUMBPALIGN=>0, MAX_SUMBPALIGN => 99999999999999,
 				MIN_AVEFRACBPMATCH =>0, MAX_AVEFRACBPMATCH => 1.01,
 				MIN_HSPS=>1, MAX_HSPS=>99999999999999,
 				NAME_TYPE => 'VERSION' ); 
	my %starts=@starts;
	my %args = (@starts,@_);
	#delete $args{''};
	foreach my $k (keys %args) {
		if (!defined $starts{$k} ) {
			print "\n$k is an invalid option for parse_query\n";
			print "  VALID OPTIONS are: ";
			foreach (keys %starts) { print "$_ ";}
			die "\nExecuation Stopped\n" ;
		}
	}
 	#print "ALLOVERLAPeq$args{'SKIP_ALLOVERLAP_ORIENT'}\n";					
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

	$line =<$FH> until $line=~ /Query= (\S+)/;
	$d{'qname'}=$1;
	if ($args{'NAME_TYPE'} eq 'VERSION') {
		$d{'qname'} = $1 if $d{'qname'} =~/\|([A-Z]+[_0-9.]+)\|/;
	}
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
 		my %sbj = ();
 		($sbj{'name'} )= $line=~/>(\S+)/;
 		if ($args{'NAME_TYPE'} eq 'VERSION') {
			$sbj{'name'} = $1 if $sbj{'name'} =~/\|([A-Z]+[_0-9.]+)\|/;
		}
 		my $tmp='';
		until ( $line =~ /Length = / ) {
			$tmp.= $line;
			$line =<$FH>;
		}
		$tmp=~s/\n/ /mg;
		$tmp=~s/ +/ /mg;
		$tmp=~s/>$sbj{'name'} *//;
		$tmp=~s/^\|//;
		$tmp=~s/ +$|^ +//m;
		$sbj{'defn'}=$tmp;
		
		
		$line=~ /= ([0-9,]+)/;
		$sbj{'len'}=$1;
		$sbj{'len'}=~s/\,//m;
		my @hsp;
		until ( $line=~ /^>/ || $line=~ /\s+Database/ ) {
			my %h=();
			$line=<$FH> until $line=~ /^\s+(Score) =/ || $line=~ /^>/ || $line=~ /\s+Database/; 
			next if $1 ne 'Score';
			( $h{'bits'}, $h{'score'} ) = $line=~ /Score = +(\S+) bits \((\d+)/;
			( $h{'expect'}) = $line =~ /Expect = (\S+)/;
			my $tmp = <$FH>;
			$line = <$FH>;
			$tmp.= $line;   #wordwrap occurs in big files with ncbi blast
			$tmp=~s/\n/ /mg;	
			($h{'bpident'},$h{'sizealign'}, $h{'%ident'}) 
					= $tmp=~/Identities = (\d+)\/(\d+) +\(([0-9.]+)/;
			($h{'bppos'},$h{'%pos'} ) = $tmp=~/Positivies = (\d+)\/\d+ +\(([0-9.]+)/;
			($h{'bpgap'}, $h{'%gap'} ) = $tmp=~/Gaps = (\d+)\/\d+ +\(([0-9.]+)/;
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
			$h{'qb'}=$qb;  $h{'qe'}=$qe;
			$h{'sb'}=$sb;	$h{'se'}=$se;
			if ($args{'ALIGN'} ) {
				$h{'alignq'}=join ("", @QA);
				$h{'alignm'}=join ("", @MA);
				$h{'aligns'}=join ("", @SA);
			}
			###additional statistics found only in this program
			$h{'bpalign'}=$h{'sizealign'} -$h{'bpgap'};
			$h{'fracbpmatch'}=$h{'bpident'}/$h{'bpalign'};
			####################################################
			###SKIP AN INDVIDUAL HSP WITH A NEXT ###############
			next if $h{'bpident'} < $args{'MIN_BPIDENT'};
			next if $h{'bpident'} > $args{'MAX_BPIDENT'};
			next if $h{'sizealign'} < $args{'MIN_SIZEALIGN'};
			next if $h{'%ident'} < $args{'MIN_%IDENT'};
			next if $h{'%ident'} > $args{'MAX_%IDENT'};
			next if $h{'bppos'} < $args{'MIN_BPPOS'};
			next if $h{'bppos'} > $args{'MAX_BPPOS'};
			next if $h{'%pos'} < $args{'MIN_%POS'};
			next if $h{'%pos'} > $args{'MAX_%POS'};
			next if $h{'bpgap'} < $args{'MIN_BPGAP'};
			next if $h{'bpgap'} > $args{'MAX_BPGAP'};
			next if $h{'%gap'} > $args{'MAX_%GAP'};
			next if $h{'%gap'} < $args{'MIN_%GAP'};		
			next if $h{'bpalign'} < $args{'MIN_BPALIGN'};
			next if $h{'bpalign'} > $args{'MAX_BPALIGN'};
			next if $h{'fracbpmatch'} < $args{'MIN_FRACBPMATCH'};
			next if $h{'fracbpmatch'} > $args{'MAX_FRACBPMATCH'};
			

			
			
			###things to bother with only if keeping query hitting self ###
			if ($args{'SKIP_SAME'} == $false) {
				#remove identical hits
				next if ( $d{'qname'} eq $sbj{'name'} ) && ($h{'qb'} == $h{'sb'}) && $args{'SKIP_IDENT_HSP'};
				#remove mirrors in for self hits
				my $skip= $false;
				foreach my $x (@hsp) {
					$skip=$true if  $h{'qb'}==$$x{'sb'} && $h{'qe'}==$$x{'se'} && $h{'sb'}==$$x{'qb'} && $h{'se'}==$$x{'qe'};
					$skip=$true if  $h{'qb'}==$$x{'se'} && $h{'qe'}==$$x{'sb'} && $h{'sb'}==$$x{'qe'} && $h{'se'}==$$x{'qb'};
				}
				next if $skip;
			}
			push @hsp, \%h;
		}
		########################################################
		#####PROCESSING BASED ON HAVING ALL HSPS AVAILABLE #####
		
		####DELETE DOUBLE OVERLAPS ##########################
		###this removes overlaps where the subject && query both overlap the other ###
		###this removes areas where overlaps can slide along repeats
		if ($args{'SKIP_OVERLAP'} && @hsp) {
			#print "DOUBLE OVERLAPS ";
			@hsp = sort { $$a{'score'} <=> $$b{'score'} } @hsp; #sort smallest to largest
			LOOPI: for (my $i=0; $i<@hsp; $i++) {
				#print "$i:$hsp[$i]{'score'}\n";
				for (my $j=$i+1;$j<@hsp; $j++) {
					if ( 		($hsp[$i]{'qb'} >= $hsp[$j]{'qb'} && $hsp[$i]{'qb'} <= $hsp[$j]{'qe'}) ||
								($hsp[$i]{'qe'} >= $hsp[$j]{'qb'} && $hsp[$i]{'qe'} <= $hsp[$j]{'qe'}) ||
								($hsp[$j]{'qb'} >= $hsp[$i]{'qb'} && $hsp[$j]{'qb'} <= $hsp[$i]{'qe'})		) {
						my ($ib,$ie, $jb, $je) = ($hsp[$i]{'sb'},$hsp[$i]{'se'},$hsp[$j]{'sb'},$hsp[$j]{'se'} );
						($ib,$ie)=($ie,$ib) if $ib > $ie;
						($jb,$je)=($je,$jb) if $jb > $je;
						if ( ($ib >= $jb && $ib <= $je) || ($ie >= $jb && $ie <= $je) || ($jb>=$ib && $jb <=$ie) ) {
							splice(@hsp, $i,1);
							redo LOOPI;
						}
					}
				}
			}
		}	
		#####REMOVE ALL OVERLAPS STARTING WITH HIGHEST SCORING PAIR ###
		#####this removes overlaps such as repeats hitting elsewhere in genome
		#####additionallly orientation must be conserved with largest hit
		if ($args{'SKIP_ALLOVERLAP_ORIENT'} && @hsp  ) {
			print "REMOVING ALL OVERLAPS: $args{'SKIP_ALLOVERLAP_ORIENT'}", scalar(@hsp)," ...\n";
			@hsp = reverse sort { $$a{$args{'SKIP_ALLOVERLAP_ORIENT'}} <=> $$b{$args{'SKIP_ALLOVERLAP_ORIENT'}} } @hsp; #sort largest to smallest
			#foreach my $h  ( @hsp ) {
			#	print "###$$h{'qb'} ($$h{'sb'})   $$h{'qe'} ($$h{'se'})  $$h{'score'}\n";
			#}
			#check if subject overlaps ##
			my $orient= substr(($hsp[0]{'se'}-$hsp[0]{'sb'}),0,1); #get - s
			#print "ORIENT: ($orient)\n";
			LOOPI: for (my $i=1; $i<@hsp; $i++) {
				####check for differing subject orientation then highest scoring pair###
				my $orient_i= substr(($hsp[$i]{'se'}-$hsp[$i]{'sb'}),0,1);
				#print "I$i($orient_i):$hsp[$i]{'qb'} ($hsp[$i]{'sb'})  $hsp[$i]{'qe'} ($hsp[$i]{'se'})    $hsp[$i]{'score'}\n";

				if ( ($orient eq "-" && $orient_i ne "-") || ($orient ne "-" && $orient_i eq "-") ) {
					splice(@hsp,$i,1);
					#print "BADORIENT";
					#my $pause=<STDIN>;
					redo LOOPI if $i < scalar(@hsp);
					last;
				
				}
				for (my $j=0;$j<$i; $j++) {
					#print "  J$j:$hsp[$j]{'qb'} ($hsp[$j]{'sb'})  $hsp[$j]{'qe'} ($hsp[$j]{'se'})    $hsp[$j]{'score'}";
					if($hsp[$i]{'qb'} >= $hsp[$j]{'qb'} && $hsp[$i]{'qb'} <= $hsp[$j]{'qe'} ) {
						#delete because query $i begins in $query $j
					} elsif($hsp[$i]{'qe'} >= $hsp[$j]{'qb'} && $hsp[$i]{'qe'} <= $hsp[$j]{'qe'} ) {
						#delete because query $i ends in $query  $j
					} elsif($hsp[$i]{'sb'} >= $hsp[$j]{'sb'} && $hsp[$i]{'sb'} <= $hsp[$j]{'se'} ) {
						#delete because subject $i begins  in subject $j
					} elsif ($hsp[$i]{'se'} >= $hsp[$j]{'sb'} && $hsp[$i]{'se'} <= $hsp[$j]{'se'} ) {
						#delete because subject $i ends  in subject $j
					} elsif ( $hsp[$i]{'qb'} > $hsp[$j]{'qe'} ) {
						#$i lies to right
						if ($orient eq '-' && ($hsp[$i]{'sb'} > $hsp[$j]{'se'}) ) {
							#bad subject on left
						} elsif ($orient ne '-' && ($hsp[$i]{'sb'} < $hsp[$j]{'se'} )) {
							#bad subject on left
						} else {
							#print "\n";
							next; #good
						}
					} elsif ( $hsp[$i]{'qe'} < $hsp[$j]{'qb'} ) {
						#$i lies to left
						if ($orient eq '-' && ($hsp[$i]{'se'} < $hsp[$j]{'sb'} ) ) {
							#bad subject lies to right
						} elsif ($orient ne '-' && ($hsp[$i]{'se'} > $hsp[$j]{'sb'} )) {
							#bad subject lies to right
						} else {
							#print "\n";
							next; #good
						}
					}
					
					#delete bad ones
					#print "DELETED";
					#my $pause=<STDIN>;
					splice(@hsp, $i,1);
					redo LOOPI if $i < scalar(@hsp);
					last;
				}
			}
			print scalar(@hsp)," left\n";
		}
		
		@{$sbj{'hsps'} } =@hsp;
		

		#####NUMBER OF HSPS AND TOTALS ########################################
		next if scalar(@hsp) < $args{'MIN_HSPS'};
		next if scalar(@hsp) > $args{'MAX_HSPS'};	

		$sbj{'sumbpalign'}=0;
		$sbj{'sumbpident'}=0;
		$sbj{'avefracbpmatch'}=0;
		foreach my $h (@hsp) { $sbj{'sumbpalign'}+=$$h{'bpalign'}; $sbj{'sumbpident'}+=$$h{'bpident'}; }
		$sbj{'avefracbpmatch'}=$sbj{'sumbpident'}/$sbj{'sumbpalign'} if $sbj{'sumbpalign'} > 0;

		next if $sbj{'sumbpalign'} < $args{'MIN_SUMBPALIGN'};
		next if $sbj{'sumbpalign'} > $args{'MAX_SUMBPALIGN'};
		next if $sbj{'avefracbpmatch'} < $args{'MIN_AVEFRACBPMATCH'};
		next if $sbj{'avefracbpmatch'} > $args{'MAX_AVEFRACBPMATCH'};
		
		###PROCESSING BASED ON THE OUTSIDE INFLUENCES ##########			
		next if  ( $d{'qname'} eq $sbj{'name'} ) && $args{'SKIP_SELF'};
		next if  $args{'SKIP_SUBJECTS'} && $args{'SKIP_SUBJECTS'}{$sbj{'name'}};
		
		push @sbjct, \%sbj
	}
	#####PROCESSING BASED ON HAVING ALL SUBJECTS AVAILABLE ####
	##########################################################
	@{ $d{'sbjct'} }=@sbjct;
 	#$line =<$FH> until $line=~/^BLAST/ || eof $FH;
  	return \%d;
}
	






######last line #####
1;
__END__

=head1 &parse_query subroutine in Blast.pm

=head2 SYNOPSIS

A call to this subroutine parses a query and returns a hash containing all of the parsed data.
Here is a short bit of code to get you started.  Pass a filehandle glob
and set any need values or filters.  Just a few of the possible filters are used in this example

 use lib '/JABPerlMod';    # directory where Blast.pm is stored
 use Blast qw(&parse_query); 
 while ( !eof(IN) ) {
    my %q = %{ &parse_query( FILEHANDLE => \*IN, 
         SKIP_SUBJECTS=> \%prev_queries,
         MIN_BPALIGN => 50, SKIP_OVERLAP=> $true, 
         SKIP_SELF => $false, SKIP_IDENT_HSP=> $true) };
    next if $q{'error'} eq 'No Blast Title Found';
    print Dumper \%q;   #use datadumper for good view of data						
    $prev_queries{$q{'qname'}}=1;
 }

=head2 DESCRIPTION

&parse_query is a subroutine that parses pairwise blastoutput.  It has many options,
but most importantly all data of interest can be parsed (including alignments).
The subroutine parses one query at a time, and will return all subjects and their 
high scoring pairs (hsps).  This program has been used extensively for nucleotide comparisons,
and has yet to be used for protein.  Additionally, it has only been used on NCBI standalone blastall
pairwise output. The parsing model was based extensively on Ian Korf's BPlite, which handles protein
and WU-BLAST, so &parse_query should be able to parse protein and WU-BLAST.  BPlite is object 
oriented which I find slightly cumbersome for quick manipulations, but may be more useful to you than
a hash structure that can be difficult to understand.

=head2 HASH STRUCTURE AND CONTENT

For each call to the subroutine, a query parsed and the parsed output is returned
in the form of a pointer to a hash.  The structure of this hash is outlined below
All values are standarded except for a few calculated 
(or value added blast) parameters denoted by ***

 
 {'<program>'}      blast program and version, e.g. BLASTN 2.0.7
 {'db'}{'name'}   name of database
 {'db'}{'seq#'}   number of sequences in database
 {'db'}{'len'}    length in bases of the database

 {'qname'}     really just the first non-space characters in title ***
                  NAME_TYPE => 'VERSION' then |VERSION| is sought
 {'qlen'}      length of query
 {'qdefn'}     definition line without name if name at beginning
 
 {'sbjct'}[0..n]  an array containing all the subject hashes
    {'name'}      name of subject (VERSION option applies) ***
    {'defn'}      defnition line without name if name at beginning
    {'len'}       length of subject sequence
    {'sumbpalign'}	sum of all hsps bpalign
    {'sumbpident'}   sum of all hsps bpident
    {'avefracbpmatch'} 	average fraction identity (sumbpident/sumbpalign)
    {'hsps'}      an array containing all of the hsp hashes
       {'bits'}       bit score for largest hsp
       {'score'}      ??? score for largest hsp
       {'expect'}     expected value based on length of query and database
       {'bpident'}    number of bases/residues identical
       {'%ident'}     % of bases/residues (bpident/sizealign)
       {'bppos'}      number of similar residues 
       {'%pos'}       % of similar residues (bppos/sizealign)
       {'bpgap'}      number of gapped bases
       {'%gap'}       % of gaps (bpgap/sizealign)
       {'sizealign'}  alignment size (bases + gaps)
       {'qb'}         beginning alignment position of query
       {'qe'}         ending alignment position of query
       {'sb'}         beginning alignment position of subject 
       {'se'}         ending alignment position of subject
       {'bpalign'}    number of bases juxtaposed (sizealign - gaps) 
       {'fracbpmatch'}bpident / bpalign   (doesn't count gaps)
       {'alignq'}     query string of alignment   (only if ALIGN => 1)
       {'alignm'}     middle matching |||||| |||| (only if ALIGN => 1)
       {'aligns'}     subject string of alignment (only if ALIGN => 1)


NOTE:A blast subject can be made of multiple hsps; however, gap blast only provides
statistics on the highest highest scoring pair.


=head2 PARSE OPTIONS

NAME_TYPE =>('')/'VERSION'    'VERSION' will parse name for version

ALIGN => (0)/1    set to 1 will return alignment strings for hsps

=head2 FILTER OPTIONS

There are going to be a lot of names that should have been deprecated before I
started programming the filters but I haven't invested the effort to add more meaningful
names yet. Below are both current and planned (***) FILTERS=>(with default values).

=head3 FILTERS for single hsp (independent of other hsps)

 ###FILTERS FOR TRADITIONAL VALUES###

 MIN_SIZEALIGN  => 0   MAX_SIZEALIGN => 99999999999
 MIN_BPIDENT    => 0   MAX_BPIDENT   => 99999999999
 MIN_%IDENT     => 0   MAX_%IDENT    => 101,
 MIN_BPPOS      => 0   MAX_BPPOS     => 99999999999
 MIN_%POS       => 0   MAX_%POS      => 101,
 MIN_BPGAP      => 0   MAX_BPGAP     => 99999999999
 MIN_%GAP       => 0   MAX_%GAP      => 101,
 #I still haven't added scoring 
 
 ###FILTERS FOR DERIVED SCORES CREATED BY THIS PROGRAM ###
 MIN_BPALIGN   (0)         MAX_BPALIGN   (9999999999)
 MIN_FRACBPMATCH (0)       MAX_FRACBPMATCH (1.01)   

=head3 FILTERS based on other hsps

SKIP_OVERLAP (0)/1       if 1 all B<double> overlaps will be removed.  Double overlaps are overlaps
where both the 2 hsps' queries overlap and the subjects overlap.  All hsps are compared 
to higher scoring ones.  If a lower scoring one overlaps it is deleted, and thus only the higher
scoring ones are keptedthe 2 hsps' queri

SKIP_ALLOVERLAP_ORIENT => (0)/1 if set to 1 all overlapping and misoriented (relative to highest hsp)
are removed.  Should be left with non-touching hsps with the same subject orientation. 

SKIP_SELF   (0)/1       if set to 1 then all hits that have the same name for query and subject
will be deleted.  Usually a record that hits itself is uninteresting.

SKIP_IDENT_HSP (0)/1     if set to 1 then all hits that have the same name for query and subject
and the same begins and ends for the alignment are deleted.  This removes a sequence hitting itself, but
allows for the seuquence to hit itself in a different location

=head3 FILTERS based on hsps totals

MIN_HSPS => (1)     MAX_HSPS => (9999999999)   The number of HSPS that must be returned 
from a subject.  MIN_HSPS => 2 is good for searching for
intron/exon structure.   

MIN_SUMBPALIGN => (0)      MAX_SUMBPALIGN => (9999999999)

MIN_AVEFRACBPMATCH =>0, MAX_AVEFRACBPMATCH => 1.01,


=head3 FILTERS for various higher level query or subject interactions

SKIP_SUBJECTS => \%HASH   A hash that contains subjects to skip.  For instance, hash 
contents of (AC00003=>1,AC00004=>1) would skip these two accession numbers if they appeared as a 
subject of a blast hit.  The value must be set to 1.  This filter is useful for removing
duplicate pairwise comparisons (A hits B and B hits A) when blasting a database against itself.
A hash is used instead of an array simply for a faster search.



=head1 AUTHOR

Jeff Bailey (jab@cwru.edu, http:)

=head1 ACKNOWLEDGEMENTS

This software was developed in the laboratory of Dr. Evan Eichler at the Department of Genetics,
Case Western Reserve University and School of Medicine
Cleveland OH 44106

=head1 COPYRIGHT

Copyright (C) 2000 Jeff Bailey. Some Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut
