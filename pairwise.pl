#!/usr/bin/env perl

#THINGS TO DO ########################################################
######################################################################

#need to add the ability to stop single stranded walks.  one sequence constant position while other keeps shifting in a tandem repeat

####THIS PROGRAM TRIMS SEARCH RESULTS IN A HEURISTIC MANNER ##########
######################################################################
#10/01/12 Tin Louie changed rm -rf to rmdir (allow another processor to finish using the subdirectory localdir/name)
#10/01/12 if -c option, then clean up /tmp/endtrim directory
#10/01/08 Tin Louie added test of the exit code of mkdir (test if another processor is working on current line)
#02=-9-06 removed the pause
#00-10-30 fixed small error in scoring and modified small score to shrink in size
#01/10/18 revamped the entire program so that it multiprocessor adept
#01/09/17 ADDED -l option to allow for loose extension when one stand reaches the end###
#01/01/13 FINALLY FIXED LACK OF ROW COLUMN LABEL
#00/09/27 FIXED OLD ERROR WITH SELF WHERE Q POS > S POS and orient S is F
#00/09/26 FIXING ERRORS IN TRACKING SYSTEM
#00/09/22 CHANGE TRACKING SO BASE NAMES LOADED#
#00/05/28 MAJOR REVAMP TO ALLOW FOR MULTIPLE PROCESSORS#

#00/05/04 FIXED THE DELETION OF THE TEMP FILES WITH RANDOM NUMBER
#00/04/04 FIXED SELF ANALYSIS AND ADDED PC FILE COMPATIBILITY (s/\r\n/\n/;)
#00/04/02 UPDATE FOR SELF ANALYSIS i.e. --nooverlap this stops any overlap
#00/04/02 ADDED RANDOM NUMBER TO QUERY AND SUBJECT AND ALIGN NAMES SO MULTIPLE COPIES CAN RUN

#!/usr/local/bin/perl
#LOAD MODULES
#use lib '/usr/people/jab/bin/JABPerlMod';
#use lib '/JABPerlMod';
#use Blast qw(&parse_query);
use Getopt::Std;
use strict 'vars';

use vars qw($true $false);
$true=1; $false=0;

use vars qw($opt_f $opt_s $opt_x  $opt_h $opt_d $opt_u $opt_o $opt_a $opt_c $opt_t $opt_v $opt_r $opt_n $opt_l);
use vars qw($puid  $pwd $tmpdir $tmp_out @files %trackbase $donedir $trackdir $localdir);
use vars qw($title $position $subseq $header $line $width_set $frac_set);
use vars qw($mismatch_score $gap_score);

use vars qw($program $pversion $pdescription $pgenerate);
$program = "$0";
$program =~ s/^.*\///;
### program stats ###
$pversion='1.011018';
$pdescription = "$program (ver:$pversion) trims the end of pairwise alignments to better define end points and to extend into repeat regions";
$pgenerate= 'jeff:dnhc eeelab: genetics:dhn';
### program stats end ###


use vars qw($default_d $default_uids $default_parm $myname);
$myname=`whoami`;
chop($myname);
$default_d="/tmp/endtrim"; 
$default_uids='endtrim:random';
$default_parm='0.3:3000:-2.3:-7';

if (!defined $ARGV[0]) {

print "USAGE
$program -t [table] -f [fasta files] [OPT ARGS]
DESCRIPTION
$pdescription
REQUIRED ARGUMENTS
-t [path] a blast parse table file is required (tab-delimited)
   *file must have first 8 columns consistant with:
    seqname1 begin1 end1 length1 seqname2 begin2 end2 length2
-f [directories] for the NON-encoded fasta files (:)
OPTIONAL ARGUMENTS
-u [buid:puid] batch unique id and program unique id
   (default $default_uids)
-o [path] output table file (default [opt -t].trim
-r [switch] start ALL over by removing all temporary directories
   *WARNING: this deletes already completed rows!!!
-c [switch] clears all tracking of other processors and crashes
   *This keeps finished rows, but will SCREW UP other
    processors still running!
-l [switch] turns off equal extension of both sequences
   *Allows extension on one sequence even if other sequence is 
    near the end and unable to add on bases for the extension process
-x [fraction:width:mismatch:gap] the trimming parameters
   *Carefully enter parameters (default $default_parm)
     fraction = the fraction of the segment width to extend
					 (recommend choices between 0.1 and 0.4)
     width    = the width in bp to anchor in alignment
                (recommend choices between 1500 and 4000)
                *the bigger the slower the program
     mismatch = the mismatch penality (negative number)
                the match penality is fixed to 1 
                (recommend choices between -2 and -3 )
     gap      = gap penality for opening gap--there is no
                penality for gap continuation
                (recommend choices between -5 and -15)
-d local directory to work in (default is /tmp/endtrim)
-s [switch] skips trimming just compiles the table
-n [switch] does not attempt to compile the table
-v [switch] more verbose output for debugging or entertainment
-a [number of lines] to see from the generated alignments
-h [switch] for help file (if it exists--currently none)
";
exit;
}

getopts('t:u:o:f:rlcnsvd:a:x:h');


if ($opt_h) {
	print "SEARCHING FOR INTERNAL DOCUMENTS...\n";
	system "perldoc $0\n";
	exit;
}

print "END TRIM **************\n";

##DEFINE THE MAIN DIRECTORY WHICH WILL INCLUDE WORKING SUBDIRECTORIES AND FINAL OUTPUTS
#$mainDir='/tmp/wgacData';

##PARSE COMMAND LINE ARGUMENTS#
$opt_t || die "Please enter with -t the path for the pairwise blast table\n";
$opt_f || die "Please enter with -f the path for the NON-encoded fasta files\n";
$opt_o || ($opt_o = $opt_t . '.trim' ) ;
$opt_l ||=0; 
$opt_d ||= $default_d;
$opt_u ||= $default_uids;
$opt_x ||= $default_parm;
$opt_l=!$opt_l ;  #reverse it so that it is normally off   no unequal is 0#
$tmp_out=$opt_o;
$tmp_out=~s/^.*\//$opt_d\//;
$pwd=`pwd`;
chomp $pwd;

#SETUP BUID AND PUID from -u#
die "Format of -u must be buid:puid\n" if $opt_u !~ /\:/;
($tmpdir,$puid) = split ":", $opt_u;
$tmpdir.="_tmp";

#SETUP THE PARMETERS FROM -x
die "Format of -x must be frac:width:mismatch:gap\n" if $opt_x !~ /\:.*\:.*\:/;
($frac_set,$width_set,$mismatch_score,$gap_score) = split ":", $opt_x;
if ($puid eq 'random') {
	$puid = "E".substr (("00000000000". int(rand(10000000000000))),-8);
}

####DEFINE THE DIRECTORIES TO BE USED #####
$trackdir="$opt_d/$tmpdir";
$donedir = "$trackdir/done";

$localdir ="$opt_d/localTrim_$tmpdir";

#if ($opt_d ) {
#	$localdir="$opt_d/$tmpdir";
#}

print "CHECKING!!! trackDir=$trackdir doneDir=$donedir localDir=$localdir opt_d=$opt_d\n";

##### RESTART COMPLETELY #########
if ($opt_r) {
	print "WE ARE STARTING ALL OVER AGAIN\n";
	if($opt_d eq $default_d){
		system "rm -rf $opt_d";
	}
	else{
		system "rm -rf $localdir";
		system "rm -rf $trackdir";
	}	
}

##### CREATE NEED DIRECTORIES #####
if($opt_d eq $default_d){
	system "mkdir -p $opt_d";
}
print "MAKING LOCAL DIRECTORIES\n";
system "mkdir -p $localdir";
system "mkdir -p $trackdir";
system "mkdir -p $donedir";

##### CLEAN UP TO REMOVED DEBRIS BUT KEEP FINISHED TRIMS ####
if ($opt_c) {
	system "rm $trackdir/processor.*";
	system "rm -rf $trackdir/row*";
	
	if ($localdir) {
		system "rm -rf $localdir/row*";
	}
}

##### READ DIRECTORIES WITHIN THE DONE DIRECTORY######
open (DONE, "find $donedir -name row* |") || die "Can't readit\n";
while (<DONE>) {
	chomp;
	my $rowname=$_;
	$rowname=~ s/^.*\///;
	$trackbase{$rowname}=$_;
}

print "ALREADY COMPLETED ",scalar( keys %trackbase) , " ROWS\n";

my $dsize=5000;


open (ACTIVE, ">$trackdir/processor.$puid\n") || die ("Can't write $trackdir/processor.$puid)!\n");
open (TABLE, $opt_t) || die "Can't open pairwise blast table: $opt_t\n";
my $header=<TABLE>;


##########################################################
###### LOOP THROUGH THE FILES ############################
$line=0;
MAIN: while ( <TABLE> ) {
	$line++;
	next if $opt_s;	
	s/\r\n/\n/;
	chomp; chomp;
	my @c=split "\t";
	
	print ACTIVE "*$line ";
	####THIS LIST IS TABDELIMITED ACCESSIONS####
	my $numb= substr ("00000000".$line,-7);
	my $name= "row$numb";
	
	if ( defined $trackbase{$name}  ) {
		print "...completed ($name).\n";
		print ACTIVE "...completed ($name).\n";
		next MAIN;
	}	
	chdir $pwd;
		
	print ACTIVE "...doing ";
	###############################################
	####READY TO RUN ##############################

	print "#####ROW $line ###($puid)###($width_set) ($frac_set)###($opt_l)#############\n";
	my $path1 = &find_file( $opt_f, '', $c[0] , '');
	my $path2 = &find_file( $opt_f, '', $c[4] , '');
	$path1="$pwd/$path1" if $path1 !~ /^\//;
	$path2="$pwd/$path2" if $path2 !~ /^\//;
	print "PATH1$path1\nPATH2$path2\n  " if $opt_v;

	my $workdir="$localdir/$name"; 
	#system "mkdir $workdir";

	my ($errorb,$errore);   #for errors defining beginning and end
	#### trim query begin and end ###
	#print "$opt_l\n";
	my ($newqb,$newsb, $newqe, $newse, $errorb, $errore)
		= &blast_trim_ends( -array => \@c ,
				-workdir=> $workdir, -fasta1=>$path1, -fasta2=>$path2, 
				-id=>$puid,-side=>'left',-width=>$width_set,
				-frac=>$frac_set, -equalextension=>	$opt_l	, -samenooverlap=>$true);
	if ($newqb eq '') {
		print "TRIMMING CRASHED... CLEANING UP AND SKIPPING!\n";
		system "rm -rf $workdir";
		next MAIN;
	}
	#### print OLD versus NEW and replace in table
	#print "OPTL:$opt_l\n";
	
	print "####OLD $c[1] ($c[5]) to $c[2] ($c[6])\n";
	print "####NEW $newqb ($newsb) to $newqe ($newse)\n";
	( $c[1],$c[5],$c[2],$c[6] )=($newqb,$newsb,$newqe,$newse);
	
	open (NEWROW, ">$donedir/$name") || die "Can't create done row ($donedir/$name)\n";
	print NEWROW join( "\t", @c),"\tB:$errorb\tE:$errore\t$name\n";
	close NEWROW;
	
	print ACTIVE "...done\n";
	system "rm -rf $workdir";
	
}
close TABLE;
close NEWTABLE;
system "mv $trackdir/processor.$puid $trackdir/$puid.log" || die "Can't move.$puid!\n";


####BAIL IF ANOTHER PROCESSOR #######
chdir "$pwd";
if ($opt_n) {
	print " This processor($puid)is quiting as forced by -n option!\n";
	exit;
}


%trackbase=();

open (FINAL, "find $donedir -name row* |") || die "Can't readit\n";
while (<FINAL>) {
	chomp;
	my $rowname=$_;
	$rowname=~ s/^.*\///;
	$trackbase{$rowname}=$_;
}
print "FINAL NUMBER COMPLETED ",scalar( keys %trackbase) , " ROWS OF $line\n";
print "COMPILING DATA INTO OUTPUT TABLE ($opt_o)!\n";

open (OUT, ">$tmp_out") || die "Can't open final output table ($tmp_out)!\n";
$header =~s/\r\n/\n/;
chomp $header;
print OUT "$header\tERRORB\tERRORE\tROW\n";


open (TABLE, $opt_t) || die "Can't open pairwise blast table: $opt_t\n";
my $header=<TABLE>;
$line=0;
while (<TABLE>) {
	s/\r\n/\n/;
	my @c=split /\t/;
	$line++;
	my $numb= substr ("00000000".$line,-7); 
	if (defined $trackbase{"row$numb"}) {
		my @header=`cat $trackbase{"row$numb"}`;
		#print $header[0];
		chomp @header;
		#print @header, "\n";
		my @r=split /\t/,$header[0];
		if ($r[0] ne $c[0] || $r[4] ne $c[4] || $r[@r-1] ne "row$numb") {
			
			print "$r[0]==$c[0]\n";
			#print "$r[1]==$c[1]\n";
			print "ROW",$r[@r-1], "\n";
			print "ERROR the following lines don't match\n";
			print "ORIGINAL", join (" ",@c)."\n";
			print "PROCESSED",join (" ",@r),"\n";
		} else {
			print OUT "$header[0]\n";
		}
	} else {
		print "ERROR row$numb file: does not exist for $line line\n";
	}
	#my $pause=<STDIN>;
}
close TABLE;
close OUT;
close ACTIVE;	

# clean up /tmp 
if ($opt_c) {

#	if($opt_d eq $default_d){
#		system "rm -rf $opt_d";
#	}
#	else{
		system "rm -rf $localdir";
		system "rm -rf $trackdir";
#	}		

}

exit;
###########################
##### END OF PROGRAM ######
###########################
	
	
#########################
#########################
#########################
######SUBROUTINES########	

sub blast_trim_ends { 		
	my %args= (-array=>'', 
				-workdir=>'', 
				-fasta1=>'',
				-fasta2=>'',
				-id=>$puid,
				-width=>3000, 
				-id=> '',
				-frac=>0.35,
				-equalextension=>0, 
				-samenooverlap=>1,
		@_);
	my $lin_count = 0;
	########################MODIFICATIONS && OTHER NOTES ##################################
	#00/04/02 UPDATE FOR SELF ANALYSIS i.e. --nooverlap this stops any overlap
	#00/04/02 ADDED RANDOM NUMBER TO QUERY AND SUBJECT AND ALIGN NAMES SO MULTIPLE COPIES CAN RUN


	#### check args and move to shorter working variable names
	die "ARG ERROR WORKDIR\n" if !$args{-workdir};
	die "ARG ERROR FASTA1\n" if !$args{-fasta1};
	die "ARG ERROR FASTA_DIR\n" if !$args{-fasta2};
	die "ARG ERROR array\n" if  !$args{-array};	
	system "mkdir -p $args{-workdir}";
	chdir $args{-workdir};
	print "WORKDIR:$args{-workdir}\n" if $opt_v;
	
	
	####################################################
	##############3 load sequences #####################
	####################################################
	my ($seq1,$head1)=&fasta_load($args{-fasta1});
	my ($seq2,$head2)=&fasta_load($args{-fasta2});
#   print"XXX",substr ($seq1,0,80),"\n";
#	print"XXX",substr ($seq2,0,80),"\n";
	
	####################################################
	
	my $ran=$args{-id};
	my $dirs=$args{'-fasta_dir'};
	my @c = @{ $args{'-array'} };
	print "$c[0] $c[1] $c[2]  $c[4]  $c[5] $c[6]\n";
	####decide if same sequence and if overlaps to be stopped#### p.s stupid not to do it
	my $stop_overlap=$false;
	$stop_overlap=$true if ($args{'-samenooverlap'} == $true && $c[0] eq $c[4]);	
   my $width=$args{'-width'};
   my $equalext=$args{-equalextension};
   print "$equalext\n";
	my ($errorb, $errore);
	#my $proportion=0.3;
	#initial calculations
	my $s_rev = $false;
   $s_rev = $true if $c[5] > $c[6];
   ###switch forward subject if self and lower positions###
	if ($c[0] eq $c[4] && $s_rev==$false && $c[1] > $c[5] && $c[1]> $c[6]) {
		 print "SWITCHERO:\n" if $opt_v;
		 ($c[5],$c[6],$c[1],$c[2])= ($c[1],$c[2],$c[5],$c[6]);
	}
	my ($newqb,$newsb,$newqe,$newse)=($c[1],$c[5],$c[2],$c[6]); 
	print "OLD:$newqb($newsb)  $newqe($newse)\n" if $opt_v;
	################################################################
	################################################# LEFT #########
	my $frac=$args{'-frac'};
	my ($same_count, $diff_count)=(0,0);
	my ($qmove_direction, $last_qmove_direction)=('','');
	my $qlast1=9999999999999;
	my $qlast2=9999999999999;
	LEFT: while (1==1) {
		#my $pause=<STDIN>;
		print "*****************************************\n" if $opt_v;
		print "**** ($newqb) \t - \t$newqe \n" if $opt_v;
		print "**** ($newsb) \t - \t$newse\n" if $opt_v;
		my $w=$width;
		my $max_w=int(($newqe-$newqb)); 
		my $max_sw=int(abs($newse-$newsb));
		$max_w=$max_sw if $max_w>$max_sw;  
		$w=$max_w if $w > $max_w; 
		print "     ID:$args{'-id'} W:$w  MAX_W:$max_w   REVERSED:$s_rev\n" if $opt_v;
		
		my ($qb,$qe)=($newqb-int($frac*$w), $newqb +$w);
		my ($sb,$se)=($newsb-int($frac*$w), $newsb +$w);
		($sb,$se)=($newsb+int($frac*$w),$newsb-$w ) if $s_rev;
		($qe,$se)=($newqe,$newse) if  $w < $width;
		 if ($qb<1) {
			my $diff= 1-$qb;
			$qb+=$diff;	#fix query
			if ($equalext) { #fix subject if equalext
				if ($s_rev) { $sb-=$diff; } else {$sb+=$diff; }
			}
		 }
		 if ($sb<1) { #reverse orientation
			my $diff=1-$sb;
			$sb+=$diff;
			$qb+=$diff if $equalext;
		 }
		 if ($sb>$c[7]) { #forward orientation
			my $diff=$sb-$c[7];  #$c[7] is length
			$sb-=$diff;
			$qb+=$diff if $equalext; #allow mismatching size  #this may need to be always on for self hits
		 }
		print "     PRE-LEFT ($frac)---$qb($sb)  $qe($se)\n";
		 if ( $newqe >= $sb && $stop_overlap && (!$s_rev)) { #overlap && stop_overlap && same_orient
		   print "PREVENTING IMPROPER OVERLAP\n" if $opt_v;
		   #I think this should still work with the $equalext option
		 	my $diff = $newqe-$sb;
		 	$sb+=$diff+1;
		 	$qb+=$diff+1;
		 }
		#print "  LEFT ($frac)   
		print "     LEFT ($frac)($newqb $newsb)---$qb($sb)  $qe($se) (SAME:$same_count) (DIFF:$diff_count)\n";
		system "rm query$ran subject$ran" if -e "query$ran";
		&fasta_save("query$ran",substr($seq1,$qb-1,$qe-$qb+1), ">query$ran");
		if ($s_rev) {
			my $rev=substr($seq2,$se-1,$sb-$se+1);
			$rev = reverse $rev;
			$rev =~ tr/atgcyrkmwsbvhdATGCYRKMWSBVHD/tacgrymkswvbdhTACGRYMKSWVBDH/; 			
			&fasta_save("subject$ran",$rev, ">subject$ran");
			
		} else {
			&fasta_save("subject$ran",substr($seq2,$sb-1,$se-$sb+1), ">subject$ran");
		}
  		print "query$ran and subject$ran and alignout$ran\n";
		system "/net/eichler/vol2/local/bin/align0 -f -30 -g -1 query$ran subject$ran  > alignout$ran";
		

		system "head -$opt_a alignout$ran	" if $opt_a;
		my ($s1,$s2)=&align_load(-alignment=>"alignout$ran");
		if ($s1 eq '') {
			print "Couldn't load alignment--bailing on this row $line\n";
			return '';

 		}
		my @p1 = align_position_array($s1);
		my @p2 = align_position_array($s2);
		my ($pb,$pe,$align_bases) = &align_pos_begin_end($s1,$s2,'left');
		print "ALIGNED BASES B:$align_bases $pb $pe\n" if $opt_v;
		#my $pause=<STDIN>;
		if ($pb==0 & $pe==0) {
			if ($width>=20000) {
				$errorb="20000 bp wouldn't align";
				last LEFT;
			}
			if ($qb==$newqb && $qe==$newqe &&  $sb==$newsb && $se ==$newse ) {
				#tryed the whole thing and it din't work
				$errorb='whole alignment bad';
				last LEFT;
			}
			print "WIDENING SEARCH FOR BEGINNING!!!\n";
			$width= $width*3;
			$width=20000 if ($width > 20000);
			next LEFT;

		}
		print "$align_bases< 0.3*($w)\n";
		if ($w < 50 ) {
			print "To small (crap alignment!\n";
			last LEFT;
		}
		if ($align_bases< 0.3*($w) && $align_bases < 1000 ) {
			if ($w==$width) {
				$width= int($width*0.8);
				$same_count=0;
				#my $pause=<STDIN>;
				next LEFT;
			}
			$errorb='minimal overlap';
		}
		my $q=$qb+$p1[$pb]-1;
		my $s=$sb+$p2[$pb]-1;
		$s=$sb-$p2[$pb]+1 if $s_rev;
		###set move###
		$last_qmove_direction=$qmove_direction;
		print "NEWBQ: $newqb ===> $q\n" if $opt_v;
		print "NEWBS: $newsb ===> $s\n" if $opt_v;
		if ($q < $newqb ) {
			$qmove_direction='<===';
		} else {
			$qmove_direction='===>';
		}
		#print "NEWB: [$q($s)] $newqe($newse) \n";
		#print "PAUSED=>"; my $pause=<STDIN>;
		if ($q==$newqb && $s==$newsb) {
			$same_count++;
			last LEFT if $frac<0.002;
			$frac=$frac*0.35;
			$diff_count=0;
		} else {
			$same_count=0;
			$diff_count++;
		}
		if ($diff_count > 5 ){
			if ($qlast2==$q) {
				print "EXACT WAFFLING 2nd to last $qlast2 same as now $q\n";
				$frac=$frac*0.8;
				$diff_count=0;
			} elsif ($qmove_direction ne $last_qmove_direction) {
				print "WAFFLING ON +/-\n";
				$frac=$frac*0.6;
				$diff_count=0;
			}
			if ($diff_count>50) {
				$frac=$frac*0.8;
				$diff_count=0;
			}
		}
		print "LAST: $last_qmove_direction NOW:$qmove_direction\n" if $opt_v;
		print "BEGIN_ALIGN_POS:$p1[$pb] (FRAC:$frac) (SAME:$same_count) (DIFF:$diff_count)\n" if $opt_v;		
		($qlast2,$qlast1)=($qlast1,$q); 
		$newqb=$q; $newsb=$s;
		#my $pause=<STDIN>;
	}
	#my $pause=<STDIN>;
	######################################################
	######################################## RIGHT #######
	$frac=$args{'-frac'};
	$width=$args{'-width'};
	$same_count=0;
	$diff_count=0;
	print "       ########\n";
	RIGHT: while(1==1) {
		#my $pause=<STDIN>;
		print "*****************************************\n" if $opt_v;
		print "**** $newqb \t - \t($newqe)\n" if $opt_v;
		print "**** $newsb \t - \t($newse)\n" if $opt_v;
		my $w=$width;
		my $max_w=int(($newqe-$newqb+1.00001)); 
		my $max_sw=int(abs($newse-$newsb)+1.00001);
		$max_w=$max_sw if $max_w>$max_sw;
		$w=$max_w if $w > $max_w;		
		print "     ID:$args{'-id'} W:$w  MAX_W:$max_w   REVERSED:$s_rev\n" if $opt_v;
		my ($qb,$qe)=($newqe-$w, $newqe +int($frac*$w));
		my ($sb,$se)=($newse-$w, $newse +int($frac*$w));
		($sb,$se)=($newse+$w,$newse- int($frac*$w)) if $s_rev;
		($qb,$sb)=($newqb,$newsb) if $w < $args{'-width'};
		if ($qe>$c[3]) { #if end #
			my $diff=$qe-$c[3];
			#print "$c[3]  $qe $diff\n";
			$qe-=$diff;
			if ($equalext) {
				if ($s_rev) { $se+=$diff;} else {$se-=$diff;}
			}
		}
		if ($se>$c[7]) { #forward orientation
			my $diff=$se-$c[7]+1;
			$se-=$diff;
			$qe-=$diff if $equalext;
		}
		if ($se<1 ) {   #reverse orientation
			my $diff = 1-$se;
			$se+=$diff;
			$qe-=$diff if $equalext;
	 	}
	 	if (   $qe >= $newsb && $stop_overlap && (!$s_rev)) { #overlap && stop_overlap && same_orient
			print "PREVENTING OVERLAP\n" if $opt_v;
			my $diff = $qe-$newsb;
			$se-=($diff+1);
			$qe-=($diff+1);
		 }

		print "     RIGHT ($frac)---($newqe $newse) $qb($sb) ($qe($se)  (SAME:$same_count) (DIFF:$diff_count)\n";
		system "rm query$ran subject$ran" if $opt_v;
		&fasta_save("query$ran",substr($seq1,$qb-1,$qe-$qb+1), ">query$ran");
		if ($s_rev) {
			my $rev=substr($seq2,$se-1,$sb-$se+1);
			$rev = reverse $rev;
			$rev =~ tr/atgcyrkmwsbvhdATGCYRKMWSBVHD/tacgrymkswvbdhTACGRYMKSWVBDH/; 			
			&fasta_save("subject$ran",$rev, ">subject$ran");
			
		} else {
			&fasta_save("subject$ran",substr($seq2,$sb-1,$se-$sb+1), ">subject$ran");
		}
		system "/net/eichler/vol2/local/bin/align0  -f -30 -g -1 query$ran subject$ran  > alignout$ran 2>/dev/null";
		system "tail -$opt_a alignout$ran	" if $opt_a;
		my ($s1,$s2)=&align_load(-alignment=>"alignout$ran");
		if ($s1 eq '') {
			print "Couldn't load alignment--bailing on this row $line\n";
			return '';
		}
		my @p1 = align_position_array($s1);
		my @p2 = align_position_array($s2);
		my ($pb,$pe,$align_bases) = &align_pos_begin_end($s1,$s2,'right');
		print "AB:$align_bases PB: $pb  PE:$pe\n" if $opt_v;
#-- modified by Lin , to stop the endless loop
		if ($pb==0 & $pe==0) {
			if ($width>=20000) {
				$errore="20000 bp wouldn't align";
				last RIGHT;
			}
			if ($qb==$newqb && $qe==$newqe &&  $sb==$newsb && $se ==$newse ) {
				#tryed the whole thing and it din't work
				$errore='whole alignment bad';
				last RIGHT;
			}
			$lin_count++;
			if($lin_count == 10){
				$lin_count=0;
				last RIGHT;
			}else{
				print "WIDENING SEARCH FOR END and $lin_count\n";
				$width*=3;
				$width=20000 if ($width > 20000);
				next RIGHT;
			}
		}
		print "$align_bases< 0.3*($w)\n";
		if ($w < 50 ) {
			print "To small (crap alignment!\n";
			last RIGHT;
		}
		if ($align_bases < 0.3*($w) && $align_bases < 1000 ) {
			if ($w==$width) {
				$width= int($width* 0.8);
				$same_count=0;
				next RIGHT;
			}
			$errore='minimal overlap';
		}
		my $q=$qb+$p1[$pe]-1;
		my $s=$sb+$p2[$pe]-1;
		$s=$sb-$p2[$pe]+1 if $s_rev;
		#for (my $x=$pe-2; $x< $pe+3;$x++) {
		#	print "$x\t",substr($s1,$x,1),'=',substr($s2,$x,1), "  $p1[$x],$p2[$x]\n";
		#}
		$last_qmove_direction=$qmove_direction;
		print "NEWB: $newqe ===> $q\n" if $opt_v ;
		if ($q < $newqe ) {
			$qmove_direction='<===';
		} else {
			$qmove_direction='===>';
		}

		#print "END_ALIGN_POS:$p1[$pe]   AB:$align_bases\n";
		#print "****NEWE: $newqb($newsb) [$q($s)]  \n";
		#print "PAUSED=>"; my $pause=<STDIN>;
		if ($q==$newqe && $s==$newse) {
			$same_count++;
			$diff_count=0;
			last RIGHT if $frac<0.002;
			$frac=$frac*0.35;
		} else {
			$same_count=0;
			$diff_count++;
		}
		if ($diff_count > 5 ){
			if ($qlast2==$q) {
				print "Exact waffling\n";
				$frac=$frac*0.5;
				$diff_count=0;
			} elsif ($qmove_direction ne $last_qmove_direction) {
				print "WAFFLING ON +/-\n";
				$frac=$frac*0.6;
				$diff_count=0;
			}
			if ($diff_count>30) {
				$frac=$frac*0.8;
				$diff_count=0;
			}
		}
		print "LAST: $last_qmove_direction NOW:$qmove_direction\n" if $opt_v;
		print "BEGIN_ALIGN_POS:$p1[$pb] (FRAC:$frac) (SAME:$same_count) (DIFF:$diff_count)\n" if $opt_v;		

		($qlast2,$qlast1)=($qlast1,$q);
		$newqe=$q; $newse=$s;
		#my $pause=<STDIN>;
	}
	unlink ("query$ran", "subject$ran", "alignout$ran");
	return $newqb, $newsb, $newqe, $newse, $errorb, $errore;
	
}


sub position_check {
	my ($qb,$qe,$sb,$se,$pc)=@_;
	my @c=@{$pc};
	my $s_rev=$false;
	$s_rev=$true if $sb > $se;

	 return $qb,$qe,$sb,$se;

}

####scoring matrux match +1 mismatch -2 gap -5####	
sub score {
	my ($b1,$b2,$gap)=@_;
	my $score=0;
	if ($b1 ne '-' && $b2 ne '-') {
		#bases
		$gap=$false;
		if ($b1 ne 'N' && $b2 ne 'N') {
			if ($b1 eq $b2) { $score=1;} else {$score= $mismatch_score;} #mismatch_socre gloabl}
		} else {
			$score= -0.0001;
		}
	} else {
		#indel
		if ($gap==$false) {
			$score= $gap_score ; #gap score is global
			$gap=$true;
		} else {
			#print "badgap\n";
			$score= -0.0001;
		}
			
	}
	return $score, $gap;
}


sub align_pos_begin_end { #0-len-1 of alignment
	my($s1,$s2, $side)=@_;
	if (length($s1) != length($s2)  ) {
		print "The align output returned did not have the same number of characters for both strands!\n";
		print "Type OK to ignore otherwise program will terminate\n";
		#my $pause=<STDIN>;
		
		#die "Sequences are different sizes!\n" if  $pause!~ /OK/i;
		if (length($s1) < length($s2) ) {
			$s1.='-' x (length($s2)-length($s1)) ;
		} else {
			$s2.='-' x (length($s1)-length($s2)) ;
		}
	}
	my $gap=$false;
	my $s;
	my ($score, $max_score, $max_position,$aligned_bases)=(0,-100000,0,0);
	#################################
	if ($side eq 'left') {
	  	############################
	  	for (my $i=length($s1)-1; $i>=0; $i--) {
		  ($s,$gap) = &score(  substr($s1,$i,1),substr($s2,$i,1),$gap );
		  $score+=$s;
		  #print "$i=>$score ($s $gap) :" if $i<800;
		 # $aligned_bases++ if $s==1; 
		  if($score>$max_score) {
			  $max_position=$i;
			  $max_score=$score;
		  }	
		}
		
	  	my $pos_begin=$max_position;
	  	print "BEGIN:$pos_begin  ($max_score)\n" if $opt_v;
	  	################################
	  	for (my $i=$pos_begin; $i<length($s1); $i++) {
		  	#print "$score " if $i<100;
		  	($s,$gap) = &score(  substr($s1,$i,1),substr($s2,$i,1),$gap );
		 	 $score+=$s;
		 	 $aligned_bases++ if $s==1; 
		  	if($score>$max_score) {
			  $max_position=$i;
			  $max_score=$score;
		 	 }
	  	}
	  	my $pos_end=$max_position;
		print "END:$pos_end  ($max_score)\n" if $opt_v;
		return $pos_begin, $pos_end, $aligned_bases;
	} else {
		#### SIDE RIGHT ###############
	  	################################
	  	###CALCULATE AND FIND MAX######
	  	for (my $i=0; $i<length($s1); $i++) {
		  	#print "$score " if $i<100;
		  	($s,$gap) = &score(  substr($s1,$i,1),substr($s2,$i,1),$gap );
		 	 $score+=$s;
		  	if($score>$max_score) {
			  $max_position=$i;
			  $max_score=$score;
		 	 }
	  	}
	  	my $pos_end=$max_position;
		#print "END:$pos_end  ($max_score)\n" if $opt_v;
		################################
		####TRACE BACK PLEASE #####
	  	for (my $i=$pos_end; $i>=0; $i--) {
		  #print "$score " if $i<100;
		  ($s,$gap) = &score(  substr($s1,$i,1),substr($s2,$i,1),$gap );
		  $score+=$s;
		 	$aligned_bases++ if $s==1;
		  if($score > $max_score) {
			  $max_position=$i;
			  $max_score=$score;
		  }	
		}
	  	my $pos_begin=$max_position;
	  	#print "BEGIN:$pos_begin  ($max_score)\n";		
	  	unlink "query";
	  	unlink "subject";
	  	unlink "alignout";
		return $pos_begin, $pos_end, $aligned_bases;	
		#################################
	}
}
sub align_position_array {
	my $s1 = shift;
	my @pos_array;
	my $position;
	for (my $i=0; $i<length($s1); $i++) {
		$position++ if substr($s1,$i,1) ne '-';
		$pos_array[$i] = $position;
		$pos_array[$i] = '-' if substr($s1,$i,1) eq '-';
	}
	return @pos_array;
}
		

sub align_load { 	
	my %args=(-alignment=>'',@_);
	my ($line, $s1,$s2);
	open (IN,$args{'-alignment'}) || return  ('','');
	$line =<IN> until $line=~/(Global) alignment score:/  || eof(IN);
	if ($1 eq "Global") {
		$line=<IN>; $line=<IN>; $line=<IN>;
		while (!eof(IN) ) {
			$line=~s/\r//;
			chomp $line;
			$line=~ / ([MBDHVRYKSWACTGN-]+)$/;
			$s1.=$1;
			$line=<IN>; $line=<IN>;
			$line=~ / ([MBDHVRYKSWACTGN-]+)$/;
			$s2.=$1;
			$line=<IN>; $line=<IN>; $line=<IN>;	 $line=<IN>;
		}
	} else {
		return ('','');
	}
	$s1=~s/ +//mg;
	$s2=~s/ +//mg;
	return $s1,$s2;
}


sub find_file {
	my ($paths, $sub_paths, $names, $extensions) = @_;
	my @paths=split ":",$paths;
	my @sub_paths=split ":",$sub_paths;
	@sub_paths=('') if !@sub_paths;
	my @names=split ":",$names;
	my @extensions=split ":",$extensions;
	@extensions=('') if ! @extensions;
	for my $path (@paths) {
	 for my $sub_path (@sub_paths) {
	 	for my $name (@names) {
	 		for my $ext (@extensions) {	
	 			my $p = $path;
	 			$p .= "/$sub_path" if $sub_path;
	 			$p .= "/$name$ext";
	 			print "TESTING $p=>" if $opt_v;
	 			if (-e $p) {print "YES\n" if $opt_v; return $p; }
	 			print "NO\n" if $opt_v;
	 		}
	 	}
	 }
	}
	print "    COULD NOT FIND FILE\n";
	return "";
}


sub fasta_load {
	open (FASTAIN, "$_[0]") || die "Can't open $_[0]\n";
	my $fasta = '';
	my $header = <FASTAIN>;
	while (<FASTAIN>) {
		s/\r\n/\n/;
		chomp;
		$fasta .= $_;
	}
	return ($fasta, $header);
	close FASTAIN;
	
}
sub fasta_save {
	my ($filename, $whole_seq, $header) = @_;
	chomp $header;
	open (FASTAOUT, ">$filename") || die "Can't create fasta ($filename)\n";
	
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
