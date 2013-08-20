#!/usr/bin/perl

###PROGRAMMING LOG (not always updated)###################################
###01-01-11 revamped program stats#
###add ability to read stdin
###add ability to output to stdout
###add ability to read its own output for -otype h########
###00-12-15 modified headers in -otype h###
###00-12-05 fix -no_query_overlap which was not removing everything#####
###00-11-28 fixed batch file processing for megablast (no repetitive database format) still need to add program record!
###00-11-26 revamp help and instructions
###00-11-26 revamp output subroutines to comply with standard nomenclature and Rsort#
###00-11-26 fixed parse command line option to allow passing of arguments#
###00-11-26 add importing of offical megablast under ba import######
###00-11-26 add additional exporting of query subject data options -ohextra -osextra#####
###00-11-09 check pp import after fixing error in preparser###
###00-11-08 fixed megablast parsing to give more accurate parsing
###00-08-20 fixed error in parser ###
###00-08-20 changed version
###00-09-19 added parsing for preparse
###00-08-13 minor modifications
###00-08-03 add improved subject exclude  and -cullsubj to &query_filter (began version2#####
###00-08-01 add blastall back to the mix to create blast_parser
###00-08-01 fixed -no_*_overlap error
###00-07-31 modified to allow for -hsort and -ssort options
###00-07-31 modified query_filter to handle fuzzy overlasp for -no_query_overlap -no_subject_overlap
###00-07-30 add documentation and do general testing (change to version 3)
###00-07-29 add documentation, alter filter subroutine, cleaning in general

###modules and pragmas#
use strict 'vars';
use Getopt::Long;
use Data::Dumper;
use Storable qw(nstore retrieve);

###declare global variables#
use vars qw ($true $false);
use vars qw (%opt @pdefaults @fdefaults %excludesubj %newdefn @files $path $ver);
use vars qw (@odefaults);
use vars qw ($excludepattern);
use vars qw ($filecount $querycount $OUT);

($true,$false)=(1,0);

###########################################################
use vars qw($program $pversion $pdescription $pgenerate);
$program = "$0";
$program =~ s/^.*\///;
### program stats ###
$pversion='42.010111';
$pdescription = "$program (ver:$pversion)  parses blast-type pairwise alignments, filters them and outputs tab-delimted tables in various formats.";
$pgenerate= 'jeff:dnhc genetics:dnh';
### program stats end ###
#print "usage: $0 -in [path] -out [path] [options]\n";
#print "DESCRIPTION\n$pdescription\n";

#html
#nullarg
#code
#description
################################################################


=head1 BLASTPARSER42

=head2 SYNOPSIS

This is a program that parses various blast-type alignments, filters them  and 
outputs them in various tabular formats.

=head2 DESCRIPTION

This is generalization of my previous NCBI blast parsers using a modular general data model that
can handle the filtering of all query/subject type nucleotide searches.  Thus, the filter options will 
be maintained from program to program.  Only the initial parsing routine to fill the query hash will
have to be changed/created for each type of program output.  There are a multitude of options which I have
strived to describe throughly within the documentation for each of the  subroutines.

Examples of blastparser uses:

C<blastparser -in file.bo -filter '-min_bpalign=\>100, -max_bpalign=\>200, -min_fracbpmatch=\>0.90'>

I<This search returns all hsps that have between 100 and 200 bases aligned with at least 90% match between them.>

=head1 THE BLOODY DETAILS

=head2 FIELDS OF THE DATA STRUCTURE BUILT OR EXPECTED BY SUBROUTINES IN THIS PROGRAM

These are the fields in the hash structure built by the parsers and expected
by the modification, filter, and output routines.  The subject and hsp field names are often
used as input for sorting and filtering routines.
Astricked subject fields are generated within the filter routine.

The following are nearly inacessible except through the code:

  {'<program>'}      blast program and version, e.g. BLASTN 2.0.7
  {'db'}{'name'}   name of database
  {'db'}{'seq#'}   number of sequences in database
  {'db'}{'len'}    length in bases of the database

  The query fields are:
  {'qname'}     really just the first non-space characters in title ***
  {'qlen'}      length of query
  {'qdefn'}     definition line without name if name at beginning
  {'sbjct'}[0..n]  an array containing all the subject hashes (NDA)

     The subject fields are:
     {'name'}      name of subject (VERSION option applies) ***
     {'defn'}      defnition line without name if name at beginning
     {'len'}       length of subject sequence
     *{'sumbpalign'}*	sum of all hsps bpalign
     *{'sumbpident'}*   sum of all hsps bpident
     *{'avefracbpmatch'}* 	average fraction identity (sumbpident/sumbpalign)
     *{'ave%ident'}*   Sum( %ident * sizealign)/ Sum(sizealign)
     {'hsps'}      an array containing all of the hsp hashes
     The hsp fields are:
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
        {'alignq'}     query string of alignment   
        {'alignm'}     middle matching string ( e.g. |||||| |||| )
        {'aligns'}     subject string of alignment 

=cut


if (!defined $ARGV[0]) {
	print "USAGE $0 -in [path]  [other options]\n";
	print "DESCRIPTION\n";
	print "$pdescription\n";
	print "PARSING\n";
	print "-in [path] to file or directory of blast/megablast/preparsed outputs\n";
	print "   STDIN uses standard input\n";
	print "-itype   blastall/megablast(ba) old megablast(mb)  preparse(pp) [default ba]\n";
	print "-parse 'OPTIONS' separated by commas[,] \n";
	print "MODIFICATIONS\n";
	print "-newscore   give a formula to be evaluated for new score\n";
	print '   default: \'$$h{fracbpmatch}**10*sqrt($$h{sumbpalign})\'',"\n";	
	print "-newdefn (path1:path2:) allows replacement of old mb's useless defnintions with better ones\n";
	print "   (the length of sequences will be hunted for in new definitions)\n";
	print "EXCLUSIONS\n";
	print "-exsubj ['pattern:string1:string2...']\n";
	print "   This pattern will be forced upon the -noprevq if selected\n";
	print "   If string1... is file then opened and  all names to exclude are read.\n";
	print "   Otherwise string1... is assumed to be a name that will be excluded\n";
	print "-noprevq [switch]  adds previous queries to -exsubj hash \n";
	print "   This switch removes redundant reporting of hits when the query becomes the subject.\n";
	print "FILTER\n";
	print "-filter 'OPTIONS' separated by commas[,]\n";
	print "OUTPUT\n";
	print "-out  [path]  path for output file (default [--in].parse)\n";
	print "-otype [letter]  (default h)\n";
	print "    (h) for hit by hit miropeat format\n";
	print "    (s) for subject by subject\n";
	print "-output OPTIONS separated by commas[,]\n";
	print "-v(erbose) show summary output when parsing\n";
	print "FULL DOCUMENTATION\n";
	print "-h(elp)\n";
	exit;
}

####GET OPTIONS###
&GetOptions(\%opt, "h","help","in:s", "itype=s","out:s" , "newdefn=s","otype:s","noprevq", "parse=s","exsubj=s", "cullsubj=s",
		"output=s", "filter=s", "osextra=s", "ohextra=s","verbose");

###LONG HELP ####
if ($opt{'h'} || $opt{'help'}) {
	system "perldoc $0\n";
	exit;
}

$opt{'in'} || die "Please designate input file with -in\n";
$opt{'itype'} ||='ba';
$opt{'itype'}='ba' if $opt{'itype'} eq 'blastall';
$opt{'itype'}='mb' if $opt{'itype'} eq 'megablast';
$opt{'otype'} ||= 'h';
$opt{'output'} ||= '';
$opt{'verbose'}||=0;
@odefaults=( );
push @odefaults, split " *[=>,]+ *",$opt{'output'};

$opt{'parse'} ||= "";
@pdefaults= ( -nametype=>'version');
push @pdefaults, split " *[=>,]+ *",$opt{'parse'};




$opt{'filter'} ||= "";
@fdefaults= ( -no_subject_self => $false);
push @fdefaults, split " *[=>,]+ *",$opt{'filter'};



$opt{'out'} || ($opt{'out'} = $opt{'in'} . '.parse');


##### GET INPUT FILE ####################
if ($opt{'in'} eq 'STDIN') {
	$path ='STDIN';
	@files=('STDIN');
} elsif (opendir (DIR, "$opt{'in'}") ) {
	@files = grep { /\w/ } readdir DIR;
	close DIR;
	$path = $opt{'in'};
} elsif (open (IN, $opt{'in'} )) {
	($path ) = $opt{'in'} =~ /(^.*)\//;
	$path ||='.';
	$opt{'in'} =~ s/^.*\///;
	@files=($opt{'in'});
} else {
	die "--in  ($opt{'in'})  is not a file and not a directory\n";
}
@files = sort @files;



#######LOAD SUBJECTS TO EXCLUDE ################
%excludesubj=();
if ( $opt{'exsubj'} ) {
	my @names;
	($excludepattern, @names) = split ":", $opt{'exsubj'};
	print "$excludepattern ", @names, "\n";
	foreach my $n (@names) {
		if (open (IN, $n) )  {
			print "FILE:($n):\n" if $opt{'verbose'};
			while (<IN>) {
				s/\r\n/\n/;
				chomp;
				/$excludepattern/;
				$excludesubj{$1}=1;
			}
		} else {
			print "NAME:($n):\n" if $opt{'verbose'};
			$n=~ /$excludepattern/;
			$excludesubj{$1}=1;
		}
	}
} 

print "SUBJECTS TO EXCLUDE: ",scalar(keys %excludesubj), "\n" if $opt{'verbose'};
my %ohash=@odefaults;
####OUTPUT HEADER LINE FORMAT ####################
$OUT='OUT';
if ($opt{'out'} eq 'STDOUT') {
	$OUT='STDOUT';
} else {
	open ( $OUT, ">$opt{'out'}") || die "Can't open out file ($opt{'out'}\n";
}
if ($opt{'otype'} eq 'h') {
	print keys %ohash, "\n" if $opt{'verbose'};
	print "$ohash{'sextra'}...$ohash{-hextra}\n" if $opt{'verbose'};
	if (lc($ohash{'-sfirst'})=~/^y/ ) {
		print $OUT "SNAME\tSB\tSE\tSLEN\tQNAME\tQB\tQE\tQLEN";
	} else {
		print $OUT "QNAME\tQB\tQE\tQLEN\tSNAME\tSB\tSE\tSLEN";
	}
	print $OUT "\tFRACBPMATCH\tBPALIGN\tSIZEALIGN\tSCORE";
	if (lc($ohash{'-sfirst'})=~/^y/ ) {print $OUT "\tSDEFN\tQDEFN";} else {print $OUT "\tQDEFN\tSDEFN";}
	if (defined $ohash{'-sextra'} ) {
		my @col=split ":",$ohash{'-sextra'};
		foreach (@col) { print $OUT "\t", uc ($_); }
	}
	if (defined $ohash{'-hextra'} ) {
		my @col=split ":",$ohash{'-hextra'};
		foreach (@col) { print $OUT  "\t", uc ($_); }
	}

	print $OUT "\n";
} elsif ($opt{'otype'} eq 's') {
	print $OUT "QNAME\tQLEN\tSDEFN\tSNAME\tSLEN\tSDEFN\t";
	print $OUT "SUMBPALIGN\tSUMBPIDENT\tAVEFRACBPMATCH\tAVEFRACBPMATCH\t";
	print $OUT "HSPS\tINDVIDUALHIT DATA=>";
	print $OUT join ":", $ohash{'-hextra'} if defined $ohash{'-hextra'};
	print $OUT "BPALIGNbp:FRACBPMATCH\%[qb-qe(sb-se)]\n";
	
} else {
	die  "Unknown output type (--otype => ($opt{'otype'})!\n";
}

#####LOAD NEW DEFNITION FILES ######################
if ($opt{'newdefn'} ) {
	my @dfiles=split ":",$opt{'newdefn'};
	my $count=0;
	foreach my $d ( @dfiles) {
		open (DEFN,$d) || die "Can't open defnition header line file ($d)!\n";
		while ( <DEFN> ) {
			chomp;
			next if !/\S/;
			my ($acc,$defn);
		   #print "$_\n";
			/>(\S+) +(.*)$/;
			($acc, $defn) = ($1,$2) ;
			#print "$acc $defn\n";
			#my $pause=<STDIN>;
			$newdefn{$acc}=$defn;
			$count++;
		}
	}
	print "TOTAL NEWDEFNS LOADED => $count\n" if $opt{'verbose'};
}
#####################################################################			
###########MAIN LOOP TO PROCESS BLAST $OUTPUT (BO) FILES #############
foreach my $file (@files) {
	$filecount++;
	my %q;
	my $FH='IN';
	if ($file ne 'STDIN') {
		open ($FH,"$path/$file") || die "Can't open input file ($path/$file)!\n";
	} else {
		$FH='STDIN';
	}
	while ( !eof($FH) ) {
		$querycount++;
		print "**** (F$filecount:Q$querycount)" if $opt{'verbose'};
		if ($opt{'itype'} eq 'mb' ) {
			print "PARSING MEGABLAST $file...\n" if $opt{'verbose'};
			%q = %{ &parse_mb_lav_single_query(-filehandle=> \*$FH, @pdefaults) };
			next if $q{'error'} eq 'No Blast Title Found';
		} elsif ($opt{'itype'} eq 'ba' ) {
			print "PARSING NCBI BLASTALL/MEGABLAST $file...\n"  if $opt{'verbose'};
			%q = %{ &parse_blastall(-filehandle=> \*$FH, @pdefaults) };
		} elsif ($opt{'itype'} eq 'pp' ) {
			print "PARSING PREPARSED BLASTALL $file...\n" if $opt{'verbose'};
			%q = %{ &parse_preparse_blastall(-filehandle=> \*$FH, @pdefaults) };
		} else {
			####add regular lav output ######
			die "-itype ($opt{'itype'}) is wrong or not yet implemented!\n";
		}
		print   "PARSE QUERY $q{'qname'} with ",scalar @{$q{'sbjct'}}," subjects...\n" if $opt{'verbose'};
		if ($opt{'newscore'} || $opt{'itype'} eq 'mb' ) {
			print "RESCORING...\n" if $opt{'verbose'};
			&rescore_the_hsps(-query => \%q );
		}
		&newdefn_replace(-query => \%q, -defn => \%newdefn );
		print "FILTERING...\n" if $opt{'verbose'};
		#print "EXCLUDEPATTERN: $excludepattern\n"; die ''
		%q=%{ &query_filter( '-query' => \%q, -exsubj=>\%excludesubj, -exsubj_pattern => $excludepattern,
					-cullsubj => $opt{'cullsubj'},  @fdefaults) };
		print "OUTPUT $q{'qname'}====[",scalar @{$q{'sbjct'}},"subjects]===>" if $opt{'verbose'};
		#print Dumper \%q;
		if ($opt{'otype'} eq 'h') {
			&print_hsps (-query=>\%q, -filehandle => \*$OUT, @odefaults);			
		} elsif ($opt{'otype'} eq 's') {
			&print_subject_per_line (-query=>\%q, -filehandle => \*$OUT, @odefaults);
		}
		if ( $opt{'noprevq'} ) {
			if ($excludepattern) {
				$q{'qname'} =~ /$excludepattern/;
				$excludesubj{$1}=1 ;
			} else {
				$excludesubj{$q{'qname'}}=1;
			}
		}

		#my $pause=<STDIN>;
		#close $OUT;
		#die "";
	}
	close $FH;
}
close $OUT;

exit ();




#####################################################################################
################### PARSING SUBROUTINES##############################################
#####################################################################################
sub parse_preparse_blastall {
 	my @starts = (-filehandle => '', #FH is glob filehandel \*IN
 		-nametype => 'version');
 		

=head2 PARSING: PREPARSE BLASTALL (-itype pp)

Parses the hsp by hsp table generated from the program preparse blastpreparse_korf.pl--which does not
currently work with megablast.

B<Options:>

-nametype => 'version', ''    version attempts to parse nbci version nomenclature from name

=cut

	my %starts=@starts;
	my %args = (@starts,@_); 
	foreach my $k (keys %args) {
		if (!defined $starts{$k} ) {
			print "\n$k is an invalid option for parse_preparse_blastall\n";
			print "  VALID OPTIONS are: ";
			foreach (keys %starts) { print "$_ ";}
			die "\nExecution Stopped\n" ;
		}
	}
	my $FH = $args{'-filehandle'}; 
	die "Expecting a glob to the filehandle e.g. FH => \\*IN\n" if $FH !~ /GLOB/;
	my %d=();
	my $line='';
	$line=<$FH> until $line=~/\t\d+\t\d+/ || $line=~/^\#/ || eof $FH;
	if ($line!~ /\t\d+\t\d+/  && $line !~/^\#/ ) {
		$d{'error'}='BAD PREPARSE DATA\n';
		return \%d;
	}
	return \%d if $line =~/\#/;
	$line=~s/\r\n/\n/;
	chomp $line;
	my @c=split "\t",$line;
	my $query=$c[0];
	$d{'qname'}= $query;
	if ($args{'-nametype'} eq 'version') {
		$d{'qname'} = $1 if $d{'qname'} =~/\|([A-Z]+[_0-9.]+)\|/;
	}
	$d{'qlen'}=$c[3];
	$d{'qdefn'}=$c[12];;
	while ($query eq $c[0]) {
		my %sbjct;
		my $subject=$c[4];
		#print "===>$subject ";
		$sbjct{'name'}=$subject;
		$sbjct{'len'}=$c[7];
		if ($args{'-nametype'} eq 'version') {
			$sbjct{'name'} = $1 if $sbjct{'name'} =~/\|([A-Z]+[_0-9.]+)\|/;
		}
		$sbjct{'defn'}=$c[13];
		while ($subject eq $c[4]) {
			
			my %h;
			$h{'qb'}=$c[1];
			$h{'qe'}=$c[2];
			$h{'sb'}=$c[5];
			$h{'se'}=$c[6];
			$h{'fracbpmatch'}=$c[8];
			$h{'bpalign'}=$c[9];
			$h{'sizealign'}=$c[10];
			$h{'score'}=$c[11];
			#11 #12 descriptions
			$h{'%ident'}=$c[14];
			$h{'bpgap'}=$c[15];
			$h{'bits'}=$c[16];
			$h{'expect'}=$c[17];
			$h{'bpident'}=$c[18];
			$h{'bppos'}=$c[19];
			###further calculations
			$h{'%gap'}=$h{'bpgap'}/$h{'sizealign'};
			$h{'%pos'}=$h{'bppos'}/$h{'sizealign'};
			push @{$sbjct{'hsps'}},\%h;
			#print ":$c[4]";
			last if eof $FH;
			$line=<$FH>;
			last if $line=~/\#/;
			$line=~s/\r\n/\n/;
			chomp $line;
			@c = split "\t", $line;
			
		}
		#print "\n";
		push @{$d{'sbjct'}},\%sbjct;
		##finish subject
		
		
		
		last if eof $FH ;
	}
		
	return \%d
}

{ my $line='';
sub parse_blastall { # -filehandle 
	#00-11-28 modified to handle megablast#
 	my @starts = (-filehandle => '', #FH is glob filehandle \*IN
   -align => 'no', 
   -nametype => 'version' ); 


=head2 PARSING: NCBI BLASTALL/MEGABLAST (-itype ba) (program default)

Parses NCBI blastall and megablast output into standard data structure.  
The only important type of data that is missing is gap number.  

B<Options:>

-align => 'no'    alignments by default are not returned if set to 'yes' 
then alignments are returned.  The alignments are not outputted in any special format.
They can be returned as part of a table.

-nametype => 'version', ''    version attempts to parse nbci version 
nomenclature from name otherwise first non-space characters are used for name.


=cut

	my %starts=@starts;
	my %args = (@starts,@_);
	#delete $args{''};
	foreach my $k (keys %args) {
		if (!defined $starts{$k} ) {
			print STDERR "\n$k is an invalid option for parse_blastall\n";
			print STDERR "  VALID OPTIONS are: ";
			foreach (keys %starts) { print STDERR "$_ ";}
			die "\nExecution Stopped\n" ;
		}
	}
	my $FH = $args{'-filehandle'}; 
 	die "Expecting a glob to the filehandle e.g. FH => \\*IN\n" if $FH !~ /GLOB/;
	 #too long type all the time
	my %d=();
	#print "$line\n";
	if ($line!~/^Query=/) {
		$line = <$FH> until $line=~/^(M?E?G?A?BLASTN? [.0-9]+)/ || eof $FH;
		if (eof $FH ) { 
			$d{'error'} = 'No Blast Title Found';
			print STDERR "ERROR NO BLAST TITLE\n";
			return \%d;
		}
		$d{'program'} = $1;
		print "P:$d{'program'}\n" if $opt{'verbose'};
		if ($d{'program'}=~/MEGA/) {
			$line=<$FH> until $line =~/^Database: +(\S+)/;
			$d{'db'}{'name'} = $1;
			$line=<$FH>;
			($d{'db'}{'seq#'}, $d{'db'}{'len'}) = $line=~ /([0-9,]+) +seq.* ([0-9,]+)/;
			$d{'db'}{'seq#'} =~ s/\,//m;
			$d{'db'}{'len'} =~ s/\,//m;	
		}
	}
	$line =<$FH> until $line=~ /Query= (\S+)/ || eof $FH;
	if (eof $FH ) { 
		$d{'error'} = 'No Hits found in megablast';
		print STDERR "ERROR NO BLAST TITLE\n";
		return \%d;
	}
	
	$d{'qname'}=$1;
	if ($args{'-nametype'} eq 'version') {
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
	$tmp=~s/\n/ /mg;
	$tmp=~s/ +/ /mg;
	$tmp=~s/^\|//;
	$tmp=~s/ +$|^ +//m;
	$tmp=~s/^Query= //;
	$d{'qdefn'}=$tmp;
	
	if ($d{'program'}=~/BLASTN/) {
		$line=<$FH> until $line =~/^Database: +(\S+)/;
		$d{'db'}{'name'} = $1;
		$line=<$FH>;
		($d{'db'}{'seq#'}, $d{'db'}{'len'}) = $line=~ /([0-9,]+) +seq.* ([0-9,]+)/;
		$d{'db'}{'seq#'} =~ s/\,//m;
		$d{'db'}{'len'} =~ s/\,//m;
	}
	#print "D:$d{'db'}{'name'}\n" if $opt{'verbose'};
	my $s=0;
	my @sbjct;
	until ( $line=~/^\s+Database/ || $line=~/^Query=/) {
		$line=<$FH> until $line=~/^(>)/ || $line=~/^\s+Database/ || $line=~/^Query=/ ;
 		next if $1 ne '>';
 		#we hit a >subject
 		my %sbj = ();
 		($sbj{'name'} )= $line=~/>(\S+)/;
 		if ($args{'-nametype'} eq 'version') {
			$sbj{'name'} = $1 if $sbj{'name'} =~/\|([A-Z]+[_0-9.]+)\|/;
		}
 		my $tmp='';
		until ( $line =~ /Length = / ) {
			$tmp.= $line;
			$line =<$FH>;
		}
		$tmp=~s/\n/ /mg;
		$tmp=~s/ +/ /mg;
		$tmp=~s/^>//;
		$tmp=~s/^$sbj{'name'} *//;
		$tmp=~s/^\|//;
		$tmp=~s/ +$|^ +//m;
		$sbj{'defn'}=$tmp;
		#print "SD:$sbj{'defn'} $line\n";
		
		$line=~ /= ([0-9,]+)/;
		$sbj{'len'}=$1;
		$sbj{'len'}=~s/\,//m;
		#print "L:$sbj{'len'}\n";
		my @hsp;
		until ( $line=~ /^>/ || $line=~ /\s+Database/ || $line=~/^Query=/ ) {
			my %h=();
			$line=<$FH> until $line=~ /^\s+(Score) =/ || $line=~ /^>/ || $line=~ /\s+Database/ || $line=~/^Query=/ ; 
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
			if ($args{'-align'} =~ /y/) {
				$h{'alignq'}=join ("", @QA);
				$h{'alignm'}=join ("", @MA);
				$h{'aligns'}=join ("", @SA);
			}
			###additional statistics found only in this program
			$h{'bpalign'}=$h{'sizealign'} -$h{'bpgap'};
			$h{'fracbpmatch'}=$h{'bpident'}/$h{'bpalign'};
			push @hsp, \%h;
		}
		@{$sbj{'hsps'} } =@hsp;		
		push @sbjct, \%sbj
	}
	@{ $d{'sbjct'} }=@sbjct;
	$line=<$FH> until $line=~/^S2:/ || $line=~/^Query=/ ;
  	return \%d;
}
} #for my $line private subroutine variable
	
###############################################################################3
#################### MEGABLAST  ROUTINES   ####################################
#############################################################################3
sub parse_mb_lav_single_query { 
####version 001108

=head2 PARSING: MB (old megablast) (-itype mb)

Parses old megablast version output.  It can not parse concatenated files--i.e.
each query must be in a separate file.  It can not currently parse -D 1 output from the 
new NCBI version of megablast.

B<Options:>

-nametype => 'version',   version attempts to parse nbci version 
nomenclature from name otherwise first non-space characters are used for name.


=cut


 	my @starts = (-filehandle => '', #FH is glob filehandle \*IN
 				-nametype => ''   #options are #version
 	); 
	my %starts=@starts;
	my %args = (@starts,@_);
	#delete $args{''};
	foreach my $k (keys %args) {
		if (!defined $starts{$k} ) {
			print STDERR "\n$k is an invalid option for parse_mb_lav_single_query\n";
			print STDERR "  VALID OPTIONS are: ";
			foreach (keys %starts) { print STDERR "$_ ";}
			die "\nExecution Stopped\n" ;
		}
	}
	my $FH = $args{'-filehandle'}; 
 	die "Expecting a glob to the filehandle e.g. FH => \\*IN\n" if $FH !~ /GLOB/;
	my %d=();
	my @sbjct=();
	my %sbj_previous=();
	my @hsp=();
	my $count=0;
	while (1==1) {
		$_=<$FH>;
		my %h=();
		if (/\#/ || eof $FH) { 	
			my @title = split("'",$_);
			my $orient = substr($title[3],0,1);
			my ($sname) = $title[1]=~/^>?(\S+)/;
			($sname) = $1 if $title[1]=~ /\|([A-Z]{1,2}\d+\.\d+)/  && $ver==$true;
			my ($qname) = $title[3]=~/[+-]?(\S+)/;
			($qname) = $1 if $title[3]=~ /\|([A-Z]{1,2}\d+\.\d+)/  && $ver==$true;	
			my ($x,$tsb,$tqb,$tse,$tqe, $score ) = split(/[\(\)\s]+/,$title[4]);
			
			#print "DATA: $sbj_previous{'name'}=>$sname $d{'qname'}=>$qname\n";
			if ($sbj_previous{'name'} ne $sname || $d{'qname'} ne $qname) {
				if (@hsp> 0) {
					@{$sbjct[$count]{'hsps'} }=@hsp;
					$sbjct[$count]{'name'}=$sbj_previous{'name'};
					$sbjct[$count]{'defn'}=$sbj_previous{'defn'};
					$count++;
					#print "PUSHING SUBJECT ",scalar(@{$d{'sbjct'}}),"\n" ;
				}
				@hsp=();
			}
			last if eof $FH;
			%h = %{ &extract_mblav() };
			if ($tqb > $tqe) {
				$h{'qb'}=$tqe; $h{'qe'}=$tqb;
				( $h{'sb'}, $h{'se'} ) = ($h{'se'},$h{'sb'} );
			}
			#print "$qname $tqb $tqe  $sname $tsb $tse\n";
			#print $h{'bpident'}/($h{'bpalign'}+$h{'gap#'});
			#print "**$qname $h{'qb'} $h{'qe'} $sname $h{'sb'} $h{se}  I:$h{'%ident'} ($h{'bpident'}/ $h{'sizealign'}) ";
			#print "ID:", $h{'bpident'}/($h{'bpalign'}+$h{'gap#'}) , "[$h{'bpident'}/($h{'bpalign'}+$h{'gap#'})] bpgap:$h{'bpgap'}\n";		
			#print "$orient $score";
			#my $pause=<STDIN>;
			#print "\n";
			#print "$qname $tqb $tqe  $sname $tsb $tse\n";
			$d{'qname'}=$qname; 
			$d{'qdefn'}=substr($title[3],1);
			$sbj_previous{'name'}=$sname;
			$sbj_previous{'defn'}=$title[1];
			push @hsp, \%h;
		}
		
	}
	###returning D####
	$d{'d'} eq '?';  #databases don't comethrough
	@{$d{'sbjct'}}=@sbjct;
	return \%d;
	
	
	die "";
	
	
	sub extract_mblav {
		#my %args=(-filehandle=>'', @_);
		#my $FH = $args{-filehandle};
		my %h=();
		my $line=<$FH>;
		#print "$line";
		die "NO a { \n" if ($line!~/^a/);
		$line=<$FH>;
		my ($nonmatch) = $line =~/(\d+)/;
		$line=<$FH>;
		($h{'sb'},$h{'qb'} ) = $line =~ /(\d+)\s+(\d+)/;
		$line=<$FH>;
		( $h{'se'},$h{'qe'} ) = $line =~ /(\d+)\s+(\d+)/;
		$line=<$FH>;
		my ($bpident,$bpgap,$bpalign,$gapnum)=(0,0,0,-1);
		my ($lqe,$lse)=(0,0);
		while ( $line=~/l (\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ ) {
			$gapnum++;
			my ($sb,$qb,$se,$qe,$perc)=($1,$2,$3,$4,$5);
			my $leng= $qe-$qb+1;
			##this is the best guess since he is doing alot of rounding##
			#print "LINE$line\n";
			if ($perc==0 || $leng<1) {
				$d{'error'}.= $line;
				print STDERR "BAD ALIGNMENT$line ($leng)\n";
				$line=<$FH>;
				next;
			}
			
			$bpalign+=$leng;
			$line=<$FH>;
			if ($lqe>0) {
				###calculate gap statistics
				my $qgap=$qb-$lqe-1;
				my $sgap=$sb-$lse-1;
				if ($qgap>0 and $sgap>0) {
					print STDERR "PRINT DOUBLE GAP ERROR\n" ;
				}
				$bpgap+=$qgap + $sgap;
			 }
			 ($lqe,$lse)=($qe,$se);
		  }
		#print "$bpgap xxxx $nonmatch\n";  	  
			my $mismatches=$nonmatch-$bpgap;
		  	  my $bpident=$bpalign-$mismatches;
		  	  my $fracbpalign=int($bpident/$bpalign*1000000)/1000000;
			  my $bptot=$bpalign+$bpgap;
			  my $percgap=int($bpgap/$bptot*10000)/100;
			  my $percident=int($bpident/$bptot*10000)/100;
			  $h{'gap#'}=$gapnum;
			  ( $h{'bpident'},$h{'sizealign'}, $h{'%ident'}, $h{'bpgap'})=($bpident,$bptot,$percident,$bpgap);
			  ( $h{'%gap'}, $h{'bpalign'}, $h{'fracbpmatch'}) = ($percgap,$bpalign,$fracbpalign);
			#print " OVERALL:$h{'score'} $h{qb}($h{sb}) to $h{qe}($h{se}) A: $h{'bpident'}/$h{'bpalign'} ($h{'fracbpmatch'})"
			#    		." I:$h{bpident}/$h{'sizealign'} ($h{'%ident'})"
			 #   		." G:$h{bpgap}/$h{'sizealign'} ($h{'%gap'}) G#:$h{'gap#'}\n";
			#my $pause=<STDIN>if $h{'fracbpmatch'} > 1;;
		return \%h;
	}
}



########################################################################
########################### MODFICATION ROUTINES #######################
########################################################################

sub newdefn_replace {

=head2 MODIFICATION replace definitions

-newdefn => path1:path2:...  The average user should never need this option.  
This option was created to allow the replacement of old mb (megablast) definitions 
with the full definitions found in the fasta files.  The length of the sequence will be hunted in the new
definitions. The length will be looked for with the pattern /L2:(\d+)/.  Definitions are load into a hash so large
number of sequences will expend large amounts of memory.

=cut


	my %args = (-query => '', -defn=> '', @_);
	my $d=$args{-defn};
	#print " OLD-queryDEFN:$args{'-query'}{'qname'} $args{'-query'}{'qdefn'}\n"; 
	if ( defined $$d{ $args{'-query'}{'qname'} } ) {
		#print "NEW -query DEFN FOUND\n";
		$args{'-query'}{'qdefn'}=$$d{ $args{'-query'}{'qname'} }
	}
	if ( $args{'-query'}{'qdefn'} =~ /L2:(\d+)/ ) { $args{'-query'}{'qlen'}=$1; }
	#print " NEW-queryDEFN:$args{'-query'}{'qname'} $args{'-query'}{'qdefn'}  $args{'-query'}{'qlen'}\n"; 
	foreach my $s ( @{ $args{'-query'}{'sbjct'} }) {
		if ($$d{$$s{'name'}} ) { $$s{'defn'}=$$d{$$s{'name'}}; }
		if ($$s{'defn'}=~/L2:(\d+)/ ) { $$s{'len'}=$1; }
			#print " NEWSUBDEFN:$$s{'name'} $$s{'defn'}  $$s{'defn'}\n"; 

	}
}


sub rescore_the_hsps {

=head2 MODIFICATION rescore hsps

-newscore => '$$h{fracbpmatch}**10*sqrt($$h{sumbpalign})'  This is default formulat.
This option must be followed by an evaluatable perl expression.
The expresion uses $$h{hsp field} for variables.  The math is limited only by perl.

=cut

	#modified 000731
	my %args = (-query => '', -formula=>'$$h{fracbpmatch} ** 5 * $$h{sizealign} ** 0.2 ', @_);
	my $scount=0;
	foreach my $s ( @{ $args{'-query'}{'sbjct'} }) {
			$scount++;
			#print "$$s{'name'}:" ;
		foreach my $h ( @{ $$s{'hsps'}   }) {
		#print "$$h{'fracbpmatch'}  and $$h{'sizealign'} \n";
		#print  eval ("$args{'-formula'}" ),"\n";
			my $newscore= eval ("$args{-formula}");
			#print "NEWSCORE $newscore\n";
			$$h{'score'}=$newscore;
			
		}
	}
	#print " ($scount)\n";
}




#############################################################################################
#####################  FILTER ROUTINE #######################################################
################################################################################################

sub query_filter { # -filehandle 
	###VERSION 1.000731###
	###00-12-19 add -all_no_query_overlap
	###00-07-31 altered no_query and no_subject_overlap to allow for some defined overlap
	###00-07-29 REVAMPED ENTIRE SUBROUTINE -lower case sub options, internal documentation,###
	


 	my @starts = (
  -query=> '', #is a parsed query input
  -exsubj=> '', #NOT VALID COMMAND LINE OPTION# 
  -exsubj_pattern => '', ##default is "blank"
  
  -no_doubleoverlap => '', #remove lower scoring hsp when both s and q overlap#
  -no_query_overlap=> '', #remove all overlapping hits in terms of query#
  -no_subject_overlap=> '', #remove all overlapping hits in terms of subject#
  -keep_single_orient => '',  #keep only orientation of best hit#
  -cullsubj =>'' ,       #remove subjects when "subject variable:pattern" (indendent of exsubj)\n";
 
  -del_self_identity=> 'no',  #when q=s remove identical hits#
  -del_self_mirror=> 'no',  #remove multiple instances of hsps for self hits#
  -no_subject_self=> 'no', #remove subjects that are same as query#

  -min_sizealign=> 0,  -max_sizealign => 9999999999, #size alignment bpalign+bpgaps#
  -min_bpident => 0,   -max_bpident =>  9999999999,  #bp matches#
  '-min_%ident' => 0, '-max_%ident' => 101,  #bpident/sizealign
  -min_bppos => 0,     -max_bppos => 9999999999,
  '-min_%pos' => 0,   '-max_%pos' =>101,
  -min_bpgap =>0,      -max_bpgap => 999999999999,
  '-min_%gap' => 0,   '-max_%gap' =>101,

  -min_bpalign=> 0,    -max_bpalign => 9999999999,    # (matches + mismatches) 
  -min_fracbpmatch => 0, -max_fracbpmatch =>1.2,   # matches / (matches + mismatches)	

  -min_sumsizealign=>0, -max_sumsizealign=>999999999999,
  -'min_ave%ident' =>0,'-max_ave%ident'=>101,      
  -min_sumbpalign=>0,   -max_sumbpalign => 99999999999999,
  -min_avefracbpmatch =>0, -max_avefracbpmatch => 1.2,
  -min_hsps=>1,  -max_hsps=>99999999999999,
  
  -all_no_query_overlap=>'',
 		); 
 		




	my %starts=@starts;
	my %args = (@starts,@_);
	foreach my $k (keys %args) {
		if (!defined $starts{$k} ) {
			print STDERR "\n$k is an invalid option for parse_query\n";
			print STDERR "  VALID OPTIONS are: ";
			foreach (sort keys %starts) { print STDERR "$_ ";}
			die "\nExecution Stopped\n" ;
		}
	}


=head2 FILTER STEP 1: REMOVAL/EXCLUSION/CULLING OF UNWANTED SUBJECTS

B<command line option>:

-exsubj ['pattern:string1:string2...']\n";
The pattern allows for less stringent matching such as getting at accession numbers even tought the actual name
is a version.  The pattern will be forced upon the -noprevq if selected\n";
The string1,string2,etc. will be treated as files first.  If the file can be opened the names will be read from the first
column of the tab-delimited file.  If the file can not be opened then the strings will be treated as names to exclude.

B<-filter option>

-cullsubj [field:pattern]  The field must be a subject field such as name.  The pattern must be
a legitimate regular expression.  If the pattern matches then the subject is excluded.  For instance 
I<name:AC> would exclude all names with accession AC.  This is independent of command line -exsubj or -exsubj_pattern.

B<subroutine info for the programmer>

-exsubj => \%HASH    This is a hash that contains excludes a subject if hash is value equal to 1. The
program user manipulates it with command line arguments -noprevq and -exsubj
For instance, hash contents of (AC00003=>1,AC00004=>1) would skip these two accession numbers if they appeared as a 
subject of a blast hit.  The value must be set to 1.  This filter is useful for removing
duplicate pairwise comparisons (A hits B and B hits A) when blasting a database against itself.
A hash is used instead of an array simply for the rapid lookup.

-exsubj_pattern => 'pattern'  The pattern allows for simpler names in the hash creation.  For instance given version
numbers for names the pattern can allow for the hash to be built from simply the accession.  So, if ver.1 is present, any other version
will be culled as well.

=cut

	#print $args{'-query'},"\n";
	my %d=%{ $args{'-query'} };
	#print $d{'qname'},"\n";
	my ($field,$pattern) = split "::" ,$args{'-cullsubj'};
	print "-cullsubj FIELD($field)PATTERN($pattern)\n" if $args{'-cullsubj'} &&  $opt{'verbose'};

	for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
		my $s= $d{'sbjct'}[$j];
		####cull subjs####################################3
		if ($args{'-cullsubj'} ) {
	
			if ( $$s{$field} =~ /$pattern/ ) {
				#print $$s{$field},"\n";

				splice( @{ $d{'sbjct'} }, $j,1);
				$j--; next;				
			}
		}
		
		####exclude subjs#################################
		##################################################
		my $exname=$$s{'name'};
		if ($args{'-exsubj_pattern'} ) {
			$exname=$1 if $$s{'name'} =~ /$args{'-exsubj_pattern'}/;
		}
		#print "EXNAME $exname  ($args{'-exsubj_pattern'})\n"; die '';
		if (${$args{'-exsubj'}}{$exname} ==1 ) {
			splice( @{ $d{'sbjct'} }, $j,1);
			$j--; next;	
		}
		###loop on hsps ###

=head2 FILTER STEP 2: single hsp (independent of other hsps)

These are the -filter options and (default) values for the single hsp filters:

 -min_sizealign  =>  (0)          -max_sizealign   => (99999999999)
 -min_bpident    =>  (0)          -max_bpident     => (99999999999)
 -min_%ident     =>  (0)          -max_%ident      => (101)
 -min_bppos      =>  (0)          -max_bppos       => (99999999999)
 -min_%pos       =>  (0)          -max_%pos        => (101)
 -min_bpgap      =>  (0)          -max_bpgap       => (99999999999)
 -min_%gap       =>  (0)          -max_%gap        => (101)
 -min_bpalign =>     (0)          -max_bpalign     => (9999999999)
 -min_fracbpmatch => (0)          -max_fracbpmatch => (1.01)   

=cut

		for( my $i =0; $i < @{$$s{'hsps'}}; $i++) {
			#print "H", scalar @{$$s{'hsps'}},"-";
			####################################################
			####REMOVE INDIVIDUAL HSPS ##########################
			#print $$s{'hsps'}[$i]{'bpident'};
			if ($$s{'hsps'}[$i]{'bpident'} < $args{'-min_bpident'}  ) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'bpident'} > $args{'-max_bpident'}  ) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'sizealign'} < $args{'-min_sizealign'} ) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'%ident'} < $args{'-min_%ident'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'%ident'} > $args{'-max_%ident'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'bppos'} < $args{'-min_bppos'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'bppos'} > $args{'-max_bppos'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'%pos'} < $args{'-min_%pos'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'%pos'} > $args{'-max_%pos'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'bpgap'} < $args{'-min_bpgap'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'bpgap'} > $args{'-max_bpgap'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'%gap'} > $args{'-max_%gap'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'%gap'} < $args{'-min_%gap'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'bpalign'} < $args{'-min_bpalign'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'bpalign'} > $args{'-max_bpalign'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'fracbpmatch'} < $args{'-min_fracbpmatch'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}
			if ($$s{'hsps'}[$i]{'fracbpmatch'} > $args{'-max_fracbpmatch'}) { splice ( @{ $$s{'hsps'} },$i,1); $i--; next;}

		}	
		if (@{ $d{'sbjct'}[$j]{'hsps'} } == 0 ) {
			splice( @{ $d{'sbjct'} }, $j,1);
			$j--; next;
		}
		
		
	}



	#####################DEAL WITH SAME HITS #####################################

=head2 FILTER STEP 3: when query and subject are the same

These are the -filter options with (defaults) for dealing with self hits:

-no_subject_self => (no)/yes       deletes subjects that have same name as query

-delete_self_mirror => (no)/yes     if set to 1 then all hits that have the same name for query and subject
and the same begins and ends for the alignment are deleted.  This removes a sequence hitting itself, but
allows for the seuquence to hit itself in a different location

-delete_self_identity => (no)/yes   removes only self identity hits 
where query postions same as subject positions

=cut

	#print "$args{'-no_subject_self'}\n";
	if (lc($args{'-no_subject_self'}) eq "yes"  ) {
		for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
			
			if ( $d{'sbjct'}[$j]{'name'} eq $d{'qname'} ) {
				#print "DELETING SELF SUBJECT ($d{'sbjct'}[$j]{'name'})\n";
				splice( @{ $d{'sbjct'} }, $j,1 );
				$j--; next;
			}
		}
	} else {
		###removes duplicate self hits ##############
		if ( lc($args{'-del_self_mirror'}) eq 'yes' ) {
			for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
				next if ( $d{'sbjct'}[$j]{'name'} ne $d{'qname'} );
				###only self_hit##
				my $s= $d{'sbjct'}[$j];
				for( my $i=0; $i < @{$$s{'hsps'}} ; $i++) {
					my $h=$$s{'hsps'}[$i];
					for(my $k=$i+1;$k<@{$$s{'hsps'}};$k++) {
						my $x=$$s{'hsps'}[$k];
						if  ( $$h{'qb'}==$$x{'sb'} && $$h{'qe'}==$$x{'se'} && $$h{'sb'}==$$x{'qb'} && $$h{'se'}==$$x{'qe'} ) {
							splice ( @{ $$s{'hsps'} },$i,1); $i--; next;
						}
						if  ($$h{'qb'}==$$x{'se'} && $$h{'qe'}==$$x{'sb'} && $$h{'sb'}==$$x{'qe'} && $$h{'se'}==$$x{'qb'} ) {
							splice ( @{ $$s{'hsps'} },$i,1); $i--; next;
						}

					}

				}
				if (@{ $d{'sbjct'}[$j]{'hsps'} } == 0 ) {splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
			}
		}
		####removes hits that are to identical positions in a self hit###
		if (lc($args{'-del_self_identity'}) eq 'yes' ) {
			for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
				next if ( $d{'sbjct'}[$j]{'name'} ne $d{'qname'} );
				###only same self_hit##
				my $s= $d{'sbjct'}[$j];
				for( my $i=0; $i < @{$$s{'hsps'}} ; $i++) {
					my $h=$$s{'hsps'}[$i];
					if  ( $$h{'qb'}==$$h{'sb'} && $$h{'qe'}==$$h{'se'} ) {
						splice ( @{ $$s{'hsps'} },$i,1); $i--; next;
					}
				}
				if (@{ $d{'sbjct'}[$j]{'hsps'} } == 0 ) {splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
			}
		}
	}
	#########################################
	########PAIRWISE COMPARISONS OF HSPS#####
   ####removes double overlaps##############
	

=head2

=cut


=head2 
FILTER STEP 4: Compare hsps within a subject

These are the -filter options for the pairwise filter options:

-no_doubleoverlap =>  [hsp field]      removes lower-valued hsps when there is overlap between 
both the query positions and the subject positions of any two hsps.  A single base pair overlap results in deletion.

-keep_single_orient => [hsp field]     removes hsps that differ in orientation to the highest scoring hsp
(this option is useful when looking for exons that should all have the same orientation)

-no_query_overlap => [hsp field]:[frac/bp]:[value]  (colon-delimited) removes lower-valued hsps when there is overlap of query positions
between any two hsps. The [hsp field] must inputted.  There is no default.
The [frac/bp]:[value] data is optional and when not specified defaults to amount of overlap being deleted. 
If [frac] then the value that follows will be a fractional amount 
of overlap in terms of the length of the higher-valued hsp from [hsp field].  Values for frac thus must be between 0 and 1.
If [bp] it will be an absolute amount of overlap in terms of number of base pairs.
Values for bp must be greater integers > 0. 

-no_subject_overlap => [hsp field]:[frac/bp]:[value] removes lower-valued hsps when there is overlap of subject positions
between any two hsps  See -no_query_overlap for [frac/bp]:[value] discussion.

=cut

	if ($args{'-no_doubleoverlap'}) {
		print "FILTERING DOUBLE OVERLAPS\n"  if $opt{'verbose'};
		for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
			next if @{ $d{'sbjct'}[$j]{'hsps'} } < 2;
			@{ $d{'sbjct'}[$j]{'hsps'} } = sort { $$a{$args{'-no_doubleoverlap'}} <=> $$b{$args{'-no_doubleoverlap'}} } @{ $d{'sbjct'}[$j]{'hsps'} };
			LOOPI: for (my $i=0; $i<@{ $d{'sbjct'}[$j]{'hsps'} }; $i++) {
				my $hi=$d{'sbjct'}[$j]{'hsps'}[$i];
				for (my $k=$i+1;$k< @{ $d{'sbjct'}[$j]{'hsps'} }; $k++) {
					my $hk=$d{'sbjct'}[$j]{'hsps'}[$k];
					#print scalar(@{ $d{'sbjct'}[$j]{'hsps'}} ), " <-hsps left  j:$j i:$i k:$k\n";
					#print "$$hi{'qb'}\t$$hi{'qe'}\t$$hi{'sb'}\t$$hi{'se'} ($$hi{'score'}) $$hi{'%ident'} $$hi{'sizealign'}\n";
					#print "$$hk{'qb'}\t$$hk{'qe'}\t$$hk{'sb'}\t$$hk{'se'} ($$hk{'score'}) $$hk{'%ident'} $$hk{'sizealign'}\n";
					#my $pause=<STDIN>;
					if ( 		($$hi{'qb'} >= $$hk{'qb'} && $$hi{'qb'} <= $$hk{'qe'}) ||
								($$hi{'qe'} >= $$hk{'qb'} && $$hi{'qe'} <= $$hk{'qe'}) ||
								($$hk{'qb'} >= $$hi{'qb'} && $$hk{'qb'} <= $$hi{'qe'})		) {
						my ($ib,$ie, $jb, $je) = ($$hi{'sb'},$$hi{'se'},$$hk{'sb'},$$hk{'se'} );
						($ib,$ie)=($ie,$ib) if $ib > $ie;
						($jb,$je)=($je,$jb) if $jb > $je;
						if ( ($ib >= $jb && $ib <= $je) || ($ie >= $jb && $ie <= $je) || ($jb>=$ib && $jb <=$ie) ) {
							splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1);
							#print "DOUBLE OVERLAP DELETING $i\n";
							redo LOOPI;
						}
					}
				}
			}
		}
	}
	#####same orientation##############
	if ($args{'-keep_single_orient'} ) { 
		for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
			#print "REMOVING ORIENTATION THAT DON'T COORESPOND TO BEST HIT: $args{'-keep_single_orientATION'}", scalar(@{ $d{'sbjct'}[$j]{'hsps'} })," ...\n";
			next if @{ $d{'sbjct'}[$j]{'hsps'} } < 2;
			@{ $d{'sbjct'}[$j]{'hsps'} } =  reverse sort { $$a{$args{'-keep_single_orientATION'}} <=> $$b{$args{'-keep_single_orientATION'}} } @{ $d{'sbjct'}[$j]{'hsps'} };

			my $orient= substr( ($d{'sbjct'}[$j]{'hsps'}[0]{'se'}-$d{'sbjct'}[$j]{'hsps'}[0]{'sb'}),0,1);
			#print "ORIENT:$orient\n";
			LOOPI: for (my $i=1; $i<@{ $d{'sbjct'}[$j]{'hsps'} }; $i++) {
				my $orient_i= substr(($d{'sbjct'}[$j]{'hsps'}[$i]{'se'}-$d{'sbjct'}[$j]{'hsps'}[$i]{'sb'}),0,1);
				#print "I:$orient_i\n";
				if ( ($orient eq "-" && $orient_i ne "-") || ($orient ne "-" && $orient_i eq "-") ) {
					splice(@{ $d{'sbjct'}[$j]{'hsps'} },$i,1);
					#print "BADORIENT";
					#my $pause=<STDIN>;
					redo LOOPI if $i < scalar(@{ $d{'sbjct'}[$j]{'hsps'} });
					last;
				}
			}
		}
	}
	if($args{'-no_query_overlap'}) {	
		my ($sortfield,$overlaptype,$amount) = split ":", $args{'-no_query_overlap'};
		#sortfield => score, begin, bpident, etc ....
		#overlaptype => frac or bp
		#amount  => either decimal or number of bases
		print "  -no_query_overlap FIELD:$sortfield TYPE:$overlaptype CUTOFF$amount\n" if $opt{'verbose'};
		for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
			next if @{ $d{'sbjct'}[$j]{'hsps'} } < 2;
			@{ $d{'sbjct'}[$j]{'hsps'} } =  reverse sort { $$a{$sortfield} <=> $$b{$sortfield} } @{ $d{'sbjct'}[$j]{'hsps'} };
			LOOPI: for (my $i=1; $i<@{ $d{'sbjct'}[$j]{'hsps'} }; $i++) {
				#my $pause=<STDIN>;
				my ($ib,$ie)= ( $d{'sbjct'}[$j]{'hsps'}[$i]{'qb'},$d{'sbjct'}[$j]{'hsps'}[$i]{'qe'} );
				#print "I:$i: $ib to $ie  $d{'sbjct'}[$j]{'hsps'}[$i]{$sortfield}\n";
				for (my $k=0;$k<$i; $k++) {
					my ($kb,$ke) = ( $d{'sbjct'}[$j]{'hsps'}[$k]{'qb'},$d{'sbjct'}[$j]{'hsps'}[$k]{'qe'} );
					#print "  K:$k: $kb to  $ke  $d{'sbjct'}[$j]{'hsps'}[$k]{$sortfield}\n";
					if($ib>=$kb && $ie <= $ke )  {
						#delete because query $i in query $k
						#print "DELETING I:$i\n";
						splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1); $i--; next LOOPI;
					}
					if ($ib<=$kb && $ie >= $ke ) {
						#delete i because query $k in query $i
						#print "DELETING I:$i\n";
						splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1); $i--; next LOOPI;
					}
					if($ib >= $kb && $ib <= $ke ) {
						#     K##################
						#                BI################
						#delete because query $i begins in $query $k
						my $bp = $ke-$ib+1;
						my $size=$ke-$kb+1;
						$size = $ie-$ib+1 if $size > ($ie-$ib+1);
						my $frac= $bp/$size;
						#print "IB IN K BP:$bp  SIZE:$size  FRAC:$frac\n";
						if ( $overlaptype eq 'frac' && $frac <= $amount) {
							#don't delete#
						} elsif ($overlaptype eq 'bp' && $bp <= $amount) {
							#don't delete#
						} else {
							#print "DELETING I:$i\n";
							splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1); $i--; next LOOPI;
						}
					}
					if($ie >= $kb && $ie <= $ke ) {
							#     					KB##################
						#              IB################IE
						#delete because query $i ends in $query  $k
						my $bp = $ie-$kb+1;
						my $size=$ke-$kb+1;
						$size = $ie-$ib+1 if $size > ($ie-$ib+1);
						my $frac= $bp/$size;
						#print "IE IN K BP:$bp  SIZE:$size  FRAC:$frac\n";
						if ( $overlaptype eq 'frac' && $frac <= $amount) {
							#don't delete#
						} elsif ($overlaptype eq 'bp' && $bp <= $amount) {
							#don't delete#
						} else {
							#print "DELETING I:$i\n";
							splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1); $i--; next LOOPI;
						}
					}
				}
			}
			#print scalar( @{ $d{'sbjct'}[$j]{'hsps'} } ), " hsps left\n";
		}	
	}
	if($args{'-no_subject_overlap'}) {		
		my ($sortfield,$overlaptype,$amount) = split ":", $args{'-no_query_overlap'};
		#sortfield => score, begin, bpident, etc ....
		#overlaptype => frac or bp
		#amount  => either decimal or number of bases
		print "\nFILTERING SUBJECT_OVERLAPS FIELD:$sortfield TYPE:$overlaptype CUTOFF$amount\n" if $opt{'verbose'};
		for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
			next if @{ $d{'sbjct'}[$j]{'hsps'} } < 2;
			@{ $d{'sbjct'}[$j]{'hsps'} } =  reverse sort { $$a{$sortfield} <=> $$b{$sortfield} } @{ $d{'sbjct'}[$j]{'hsps'} };
			LOOPI: for (my $i=1; $i<@{ $d{'sbjct'}[$j]{'hsps'} }; $i++) {
				my ($ib,$ie)= ( $d{'sbjct'}[$j]{'hsps'}[$i]{'sb'},$d{'sbjct'}[$j]{'hsps'}[$i]{'se'} );
				($ib,$ie)=($ie,$ib) if $ib>$ie;
				#print "I:$i: $ib to $ie  $d{'sbjct'}[$j]{'hsps'}[$i]{$args{'-no_subject_overlap'}}\n";
				for (my $k=0;$k<$i; $k++) {
					my ($kb,$ke) = ( $d{'sbjct'}[$j]{'hsps'}[$k]{'sb'},$d{'sbjct'}[$j]{'hsps'}[$k]{'se'} );
					($kb,$ke)=($ke,$kb) if $kb>$ke;
					#print "  K:$k: $kb to  $ke  $d{'sbjct'}[$j]{'hsps'}[$k]{$args{'-no_subject_overlap'}}\n";
					if($ib>=$kb && $ie <= $ke )  {
						#delete because query $i in query $k
						#print "DELETING I:$i\n";
						splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1); $i--; next LOOPI;
					}
					if ($ib<=$kb && $ie >= $ke ) {
						#delete i because query $k in query $i
						#print "DELETING I:$i\n";
						splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1); $i--; next LOOPI;
					}
					if($ib >= $kb && $ib <= $ke ) {
						#     K##################
						#                BI################
						#delete because query $i begins in $query $k
						my $bp = $ke-$ib+1;
						my $size=$ke-$kb+1;
						$size = $ie-$ib+1 if $size > ($ie-$ib+1);
						my $frac= $bp/$size;
						#print "IB IN K BP:$bp  SIZE:$size  FRAC:$frac\n";
						if ( $overlaptype eq 'frac' && $frac <= $amount) {
							#don't delete#
						} elsif ($overlaptype eq 'bp' && $bp <= $amount) {
							#don't delete#
						} else {
							#print "DELETING I:$i\n";
							splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1); $i--; next LOOPI;
						}
					}
					if($ie >= $kb && $ie <= $ke ) {
							#     					KB##################
						#              IB################IE
						#delete because query $i ends in $query  $k
						my $bp = $ie-$kb+1;
						my $size=$ke-$kb+1;
						$size = $ie-$ib+1 if $size > ($ie-$ib+1);
						my $frac= $bp/$size;
						#print "IE IN K BP:$bp  SIZE:$size  FRAC:$frac\n";
						if ( $overlaptype eq 'frac' && $frac <= $amount) {
							#don't delete#
						} elsif ($overlaptype eq 'bp' && $bp <= $amount) {
							#don't delete#
						} else {
							#print "DELETING I:$i\n";
							splice(@{ $d{'sbjct'}[$j]{'hsps'} }, $i,1); $i--; next LOOPI;
						}
					}

				}
			}
			#print scalar( @{ $d{'sbjct'}[$j]{'hsps'} } ), " hsps left\n";
		}	
	}
	#################################################################
	#################################################################
	####PARSING BASED ON HSPS REMAINING #############################


=head2 FILTER STEP 5: subject hsp totals (all hsps)

These are the -filter options with (defaults) for the sum of subject hsps.

 -min_sumsizealign =>     (0)   -max_sumsizealign =>   (99999999999999)
 -min_ave%ident =>        (0)   -max_ave%ident =>      (101)
 -min_sumbpalign =>       (0)   -max_sumbpalign =>     (99999999999999)
 -min_avefracbpmatch =>   (0)   -max_avefracbpmatch => (1.2)
 -min_hsps =>             (1)   -max_hsps=>            (99999999999999)

=cut
	
	for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
		if ( @{ $d{'sbjct'}[$j]{'hsps'} } < $args{'-min_hsps'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
		if ( @{ $d{'sbjct'}[$j]{'hsps'} } > $args{'-max_hsps'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
		$d{'sbjct'}[$j]{'sumbpalign'}=0; 
		$d{'sbjct'}[$j]{'sumbpident'}=0;
		
		$d{'sbjct'}[$j]{'avefracbpmatch'}=0;
		$d{'sbjct'}[$j]{'ave%ident'}=0;
		$d{'sbjct'}[$j]{'sumsizealign'}=0;
		
		foreach my $h (@{ $d{'sbjct'}[$j]{'hsps'} } ) { 
			$d{'sbjct'}[$j]{'ave%ident'}+= $$h{'sizealign'}*$$h{'%ident'};
			$d{'sbjct'}[$j]{'sumsizealign'}+= $$h{'sizealign'};		
		
			$d{'sbjct'}[$j]{'sumbpalign'}+= $$h{'bpalign'}; 
			$d{'sbjct'}[$j]{'sumbpident'}+= $$h{'bpident'};
		}
		$d{'sbjct'}[$j]{'ave%ident'}=$d{'sbjct'}[$j]{'ave%ident'}/$d{'sbjct'}[$j]{'sumsizealign'} if $d{'sbjct'}[$j]{'sumsizealign'} >0;
		$d{'sbjct'}[$j]{'avefracbpmatch'}=$d{'sbjct'}[$j]{'sumbpident'}/$d{'sbjct'}[$j]{'sumbpalign'} if $d{'sbjct'}[$j]{'sumbpalign'} > 0;
		if ( $d{'sbjct'}[$j]{'sumsizealign'} < $args{'-min_sumsizealign'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
		if ( $d{'sbjct'}[$j]{'sumsizealign'} > $args{'-max_sumsizealign'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
		if ( $d{'sbjct'}[$j]{'ave%ident'} < $args{'-min_ave%ident'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
		if ( $d{'sbjct'}[$j]{'ave%ident'} > $args{'-max_ave%ident'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}		
		
		if ( $d{'sbjct'}[$j]{'sumbpalign'} < $args{'-min_sumbpalign'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
		if ( $d{'sbjct'}[$j]{'sumbpalign'} > $args{'-max_sumbpalign'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
		if ( $d{'sbjct'}[$j]{'avefracbpmatch'} < $args{'-min_avefracbpmatch'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
		if ( $d{'sbjct'}[$j]{'avefracbpmatch'} > $args{'-max_avefracbpmatch'} ) { splice( @{ $d{'sbjct'} }, $j,1);$j--; next;}
	}

=head2 FILTER STEP 6: pairwise comparisons between all hsps

=cut
	if($args{'-all_no_query_overlap'}) {	
		my ($scorefield,$overlaptype,$amount) = split ":", $args{'-all_no_query_overlap'};
		#scorefield => score, begin, bpident, etc ....
		#overlaptype => frac or bp
		#amount  => either decimal or number of bases
		print "  -all_no_query_overlap FIELD:$scorefield TYPE:$overlaptype CUTOFF$amount\n" if $opt{'verbose'};
		my @h=();
		for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
			@{ $d{'sbjct'}[$j]{'hsps'} } =  sort { $$a{'qb'} <=> $$b{'qb'} || $$a{'qe'} <=> $$b{'qe'} } @{ $d{'sbjct'}[$j]{'hsps'} };
			$h[$j]=0;
		}
		
		#my $pause=<STDIN>;
		WLOOP: while (1==1) {
			#print "ARRAY#",scalar(@{ $d{'sbjct'} }),"\n";
			my $more=$false;
			#####find the the query with lowest qb####
			my $low_s='NA';
			print "ARRAY#",scalar(@{ $d{'sbjct'} }),"\n";
			for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
				print "J$j(",scalar @{ $d{'sbjct'}[$j]{'hsps'} }, ")";
				if ($h[$j] >= @{ $d{'sbjct'}[$j]{'hsps'} }) {
					print "\n";
					next;
				}
				$low_s=$j if $low_s eq 'NA';
				$more=$true;
				
				print "==>$h[$j]=>$d{'sbjct'}[$j]{'hsps'}[$h[$j]]{'qb'}-$d{'sbjct'}[$j]{'hsps'}[$h[$j]]{'qe'} $d{'sbjct'}[$j]{'hsps'}[$h[$j]]{'fracbpmatch'}\n";
				$low_s=$j if $d{'sbjct'}[$j]{'hsps'}[$h[$j]]{'qb'} < $d{'sbjct'}[$low_s]{'hsps'}[$h[$low_s]]{'qb'} ;
			}
			last WLOOP if $more==$false;
			my $lh=$d{'sbjct'}[$low_s]{'hsps'}[$h[$low_s]];
			print "LOWS:$low_s=>$h[$low_s]=>$$lh{'qb'} $$lh{'qe'} $$lh{'frabpmatch'}\n";
			my $no_overlap=$true;
			#my $pause=<STDIN>;
			FLOOP: for (my $j=0; $j< @{ $d{'sbjct'} }; $j++ ) {
				next if $h[$j] >= @{ $d{'sbjct'}[$j]{'hsps'} };
				next if $low_s==$j;
				if ($d{'sbjct'}[$j]{'hsps'}[$h[$j]]{'qb'} <= $$lh{'qe'}) {
					####we have an overlap to resovlve####
					my $jh=$d{'sbjct'}[$j]{'hsps'}[$h[$j]];
					my $bpoverlap = $$lh{'qe'}-$$jh{'qb'}+1;
					my $smallbpalign = $$lh{'bpalign'};
					$smallbpalign = $$jh{'bpalign'} if $$jh{'bpalign'} < $$lh{'bpalign'};
					my $fracbpoverlap=$bpoverlap/$smallbpalign;
					print "OVERLAP=>LOW $low_s:$h[$low_s] ($$lh{'qb'}-$$lh{'qe'} $$lh{$scorefield}) $j:$h[$j]($$jh{'qb'}-$$jh{'qe'} $$jh{$scorefield}) $bpoverlap \n";
					if ( $overlaptype eq 'frac' && $fracbpoverlap <= $amount) {
						#don't delete#
						
					} elsif ($overlaptype eq 'bp' && $bpoverlap <= $amount) {
						#don't delete#
					} else {
						#choose which to delete#
						$no_overlap=$false;
						
						print "J$$jh{$scorefield} <=> L$$lh{$scorefield}\n";
						my $todel= ($$jh{$scorefield} <=> $$lh{$scorefield});
						print "TODEL_SIGN:$todel\n";
						if ($todel>0) {
							print "DEL LOW_S\n";
							$todel=$low_s;
						} elsif( $todel <0 ) {
							print "DEL J\n";
							$todel=$j;
						} else {
							$todel=$j;
							$todel=$low_s if $h[$low_s]<=$h[$j];
						}	
						print "   Delete DEL$todel:H$h[$todel] N",scalar @{$d{'sbjct'}[$j]{'hsps'}},"\n";
						splice(@{ $d{'sbjct'}[$todel]{'hsps'} }, $h[$todel],1); 
						if (@{ $d{'sbjct'}[$todel]{'hsps'} }==0) {
							splice(@{ $d{'sbjct'} },$todel,1);
							splice(@h,$todel,1);
						}
						next WLOOP;

					}
				}
			}
			print "NO OVERLAP INCREMENTING $low_s\n";
			$h[$low_s]++ if $no_overlap==$true;
			#my $pause=<STDIN>;
		}
		print "\n";
	}

 
	#print "RETURNING D\n\n";
	return \%d;


}	


sub print_hsps {
	#modified 001126


=head2 OUTPUT: Hsp by Hsp (-otype h)

This -otype prints a summary line for each hsp (pairwise comparison).  This output is compatable with miropeats output of the standard:
name1:begin1:end1:len1:name2:begin2:end2:len2.  

B<Options:>

-sfirst=>(no)/yes  Default is no. When yes the subject will be switched in position with the query.  The subject will always have a positive
orientation while the query may now be reversed.

-ssort=>name  The subject field used to sort the subjects.  As this occurs after filtering, subject fields that
are generated during filtering such as avefracbpmatch may be used.  The default is query name.

-hsort=>qb  The hsp field used to sort the hsps (default qb).  Any numerical hsp field may be used.

-hextra=>''  Additional hsp field(s) delimited by colons(:) to be outputted.

-sextra=>''  Additional subject field(s) delimited by colons(:) to be outputted.

Note that -sextra does not exist because currently all valid subject fields are already in output.


=cut


	#########DEFAULTS####################
	my @starts = (-filehandle => \*STDOUT, -query => '', 
   -sfirst=>'no', #'yes' puts suject info first 
   -ssort=>'name', #subject field to sort
   -hsort=>'qb',   #hsp field to sort 
   -hextra=>'',   #additonal hsp fields to add to output
   -sextra=>''  #additonal subject fields to add to output
  	);
   my %starts=@starts;
	my %args = (@starts,@_); 
	foreach my $k (keys %args) {
		if (!defined $starts{$k} ) {
			print STDERR "\n$k is an invalid option for -output for output type h!\n";
			print STDERR "  VALID OPTIONS are: ";
			foreach (keys %starts) { print STDERR "$_ ";}
			die "\nExecution Stopped\n" ;
		}
	}
	my $fh= $args{'-filehandle'};
	my $hsort=$args{'-hsort'};
	my $ssort=$args{'-ssort'};
	$args{'-sfirst'}=lc $args{'-sfirst'};
	#print "SFIRST:$args{'-sfirst'}\n";
	my @osextra=();
	@osextra=split ":", $args{'-sextra'} if $args{'-sextra'};
	my @ohextra=();
	@ohextra=split ":", $args{'-hextra'} if $args{'-hextra'};
	
	my $scount=0;
	my $hcount=0;
	my $hrev=$false;
	$hrev=$true if $hsort =~ s/^R// ;
	my $srev=$false;
	$srev=$true if $ssort =~s/^R//;
	
	my @sorts=sort {$$a{$ssort} <=> $$b{$ssort} } ( @{ $args{'-query'}{'sbjct'} });
	@sorts=reverse @sorts if $srev;
	foreach my $s ( @sorts) {
			$scount++;
			print "$$s{'name'} "  if $opt{'verbose'} ;
		my @sorth= sort { $$a{$hsort} <=> $$b{$hsort} }  (@{ $$s{'hsps'}   });
		@sorth=reverse @sorth if $hrev;
		foreach my $h ( @sorth ) {
			$hcount++;
			if ($args{'-sfirst'} =~/^y/) {
				if ( $$h{sb} < $$h{se} ) {
					print $fh "$$s{name}\t$$h{sb}\t$$h{se}\t$$s{len}";
					print $fh "\t$args{'-query'}{qname}\t$$h{qb}\t$$h{qe}\t$args{'-query'}{'qlen'}";
				} else {
					print $fh "$$s{name}\t$$h{se}\t$$h{sb}\t$$s{len}";
					print $fh "\t$args{'-query'}{qname}\t$$h{qe}\t$$h{qb}\t$args{'-query'}{'qlen'}";				
				}
			} else {
				print $fh "$args{'-query'}{qname}\t$$h{qb}\t$$h{qe}\t$args{'-query'}{'qlen'}";
				print $fh "\t$$s{name}\t$$h{sb}\t$$h{se}\t$$s{len}";
			} 
			#print "($$h{'fracbpmatch'})\n";
			print $fh "\t$$h{'fracbpmatch'}\t$$h{'bpalign'}\t$$h{'sizealign'}";
			print $fh "\t$$h{'score'}";
			if ($args{'-sfirst'}=~/^y/) {
				print $fh "\t$$s{'defn'}\t$args{'-query'}{qdefn}";
			} else {
				print $fh "\t$args{'-query'}{qdefn}\t$$s{'defn'}";
			}
			foreach (@osextra) {print $fh "\t$$s{$_}";}
			foreach (@ohextra) {print $fh "\t$$h{$_}";}
			print $fh "\n";
			
		}
	}
	print "\nTOTAL (SUB $scount:HSPS $hcount)\n" if $opt{'verbose'};
}
sub print_subject_per_line {

=head2 $OUTPUT: Subject Summary (-otype s)

This -otype prints a summary of a single subject per line.

B<Options:>

-ssort=>'name'  The subject field used to sort the subjects.  As this occurs after filtering, subject fields that
are generated during filtering such as avefracbpmatch may be used.  The default is query name.

-hsort=>'sb',   The hsp field use to sort the hsps (default sb).  Any numerical hsp field may be used.

-hextra=>''  Additional hsp field(s) delimited by colons(:) to be outputted.

Note that -sextra does not exist because currently all valid subject fields are already in output.


=cut



	my @starts = (-filehandle => \*STDOUT, -query => '', 
   -ssort=>'name', #subject field to sort 
   -hsort=>'sb',   #hsp field to sort   
   -hextra=>'', #additonal fields to add
   );
	my %starts=@starts;
	my %args = (@starts,@_); 
	foreach my $k (keys %args) {
		if (!defined $starts{$k} ) {
			print STDERR "\n$k is an invalid option in -output for output type s!\n";
			print STDERR "  VALID OPTIONS are: ";
			foreach (keys %starts) { print  STDERR "$_ ";}
			die "\nExecution Stopped\n" ;
		}
	}
	my $fh= $args{'-filehandle'};
	my $hsort=$args{'-hsort'};
	my $ssort=$args{'-ssort'};
	my @ohextra=();
	@ohextra=split ":", $args{'-hextra'} if $args{'-hextra'};
	my $scount=0;
	my $hcount=0;
	my $hrev=$false;
	$hrev=$true if $hsort =~ s/^R// ;
	my $srev=$false;
	$srev=$true if $ssort =~s/^R//;
	
	my @sorts=sort {$$a{$ssort} <=> $$b{$ssort} } ( @{ $args{'-query'}{'sbjct'} });
	@sorts=reverse @sorts if $srev;
	foreach my $s ( @sorts) {
		$scount++;
		print "$$s{'name'}:"  if $opt{'verbose'};
		print $fh "$args{'-query'}{qname}\t$args{'-query'}{'qlen'}\t$args{'-query'}{'qdefn'}\t";
		print $fh "$$s{name}\t$$s{len}\t$$s{'defn'}\t";
		print $fh "$$s{'sumbpalign'}\t$$s{'sumbpident'}\t$$s{'avefracbpmatch'}\t$$s{'avefracbpmatch'}\t";
		print $fh scalar (@{ $$s{'hsps'}   });
		
		my @sorth= sort { $$a{$hsort} <=> $$b{$hsort} }  (@{ $$s{'hsps'}   });
		@sorth=reverse @sorth if $hrev;
		foreach my $h ( @sorth ) {
			$hcount++;
			print $fh "\t";
			foreach (@ohextra) { print $fh "$$h{$_}:"; }
			print $fh "$$h{'bpalign'}bp:",int(10000*$$h{'fracbpmatch'})/100, "\%[$$h{qb}-$$h{qe}($$h{sb}-$$h{se})]";
			
		}
		print $fh "\n";

	}
	print "\nTOTAL (SUB $scount:HSPS $hcount)\n" if $opt{'verbose'};
}




=head1 AUTHOR

Jeff Bailey (jab@cwru.edu, http:)

=head1 ACKNOWLEDGEMENTS

This software was developed in the laboratory of:
 Dr. Evan Eichler 
 Department of Genetics,
 Case Western Reserve University and School of Medicine
 Cleveland OH 44106

=head1 COPYRIGHT

Copyright (C) 2000 Jeff Bailey. Extremely Large Monetary Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut
