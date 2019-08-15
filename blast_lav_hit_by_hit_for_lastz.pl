#!/usr/bin/perl
#THIS PROGRAM PARSES BLASTOUTPUT FROM WEBB MILLER'S LAV OUTPUT IN A HIT BY HIT BASIS
#00/04/04 reactivated "what to do if skip_same" is false (this checks for mirror images in self hits#
#00/04/02 error fixed for empty strings [next if $bpalign <=0; ] added to bad alignments coming from lav_self_overlap_breaker2.pl
use strict 'vars';

use Getopt::Long;
use Data::Dumper;
use FindBin;
# changed DG
use lib "./JABPerlMod";
#use lib "$FindBin::Bin/JABPerlMod"; # subdirectory JABPerlMod should contain Blast.pm
#use lib "$FindBin::Bin/JABPerlMod"; # subdirectory JABPerlMod should contain Blast.pm
# end changed DG
use Blast qw(&parse_query);

use vars qw($true $false);
$true=1; $false=0;


use vars qw(%opt);
use vars qw(%prev_queries @files $path @defaults);


if ($ARGV[0] eq '') {
	print "blast_lav_hit_by_hit.pl
  --in [file/directory] containing blast records to parse
  --out [filename] for output table (default [in].parse)
  --(no)header [switch] for output table (default header)
  --(no)previousq [switch] remove previous queries from being subjects (default nopreviousq)
  --options
  ";
	die "";
}
####GET OPTIONS
&GetOptions(\%opt, "in=s", "out:s", "previousq!", "header!","options:s");
$opt{'in'} || die "Please designate input file with -in\n";
$opt{'options'} || ( $opt{'options'} = "" ) 	;
$opt{'options'} = ', '.$opt{'options'} if $opt{'options' ne ''};
$opt{'out'} || ($opt{'out'} = $opt{'in'} . '.parse');
 ($opt{'header'} eq '0') || ($opt{'header'}=1);
@defaults= (SKIP_OVERLAP=> $true, SKIP_SELF => $false, SKIP_IDENT_HSP=> $true);
	#		'MAX_%GAP'=> 80, MIN_BPALIGN => 10, MIN_FRACBPMATCH=>0.80,
	#		,SKIP_OVERLAP=> $false, SKIP_SELF => $false, SKIP_IDENT_HSP=> $false) };
push @defaults, split " *[=>,]+ *",$opt{'options'};

%prev_queries=();
if (opendir (DIR, $opt{'in'}) ) {
	@files = grep { /[a-zA-Z0-9]/ } readdir DIR;
	@files = sort @files;
	close DIR;
	$path = $opt{'in'};
} elsif (open (IN, $opt{'in'})) {
	($path ) = $opt{'in'} =~ /(^.*)\//;
	$opt{'in'} =~ s/^.\///;
	@files=($opt{'in'});
} else {
	die "--in  ($opt{'in'})  not a file and not a directory\n";
}

# added error handling if can't open the file.  
# David Gordon, Aug 25, 2016

open (OUT, ">$opt{'out'}") || die "couldn't open " . $opt{'out'} . " for write\n";
if ($opt{'header'}) {
	print OUT "NAME1\tBEGIN1\tEND\tLEN1\tNAME2\tBEGIN2\tEND2\tLEN2";
	print OUT "\tFRAC_BP_IDENT\tBP_ALIGN\tSIZE_ALIGN\tSCORE";
	print OUT "\tDEFN1\tDEFN2";
	print OUT "\n";
}
foreach my $f (@files) {
	print "$path/$f\n";
	open (IN, "$path/$f") || die "Can't open file ($f)! \n";
	my %q = %{ &parse_lav( FILEHANDLE => \*IN , SKIP_SUBJECTS=> \%prev_queries,
					@defaults)};

	next if $q{'error'} eq 'No Blast Title Found';
	print "$q{'qname'}=====>";
	#print Dumper \%q;
	&parse_query_format_pairwise (QUERY=>\%q, FILEHANDLE => \*OUT);				
						
	$prev_queries{$q{'qname'}}=1;
}

sub parse_lav{ # FILEHANDLE 
 	my %args = (FILEHANDLE => '', #FH is glob filehandle \*IN
 				ALIGN => $false, SKIP_OVERLAP => $false, 
 				SKIP_ALLOVERLAP_ORIENT => $false,
 				SKIP_IDENT_HSP=> $false, SKIP_SELF=> $false,
 				'MIN_SIZEALIGN'=> 0, 
 				'MIN_BPIDENT' => 0,
 				'MIN_%IDENT' => 0, 'MAX_%IDENT' => 101,
 				'MIN_BPOS' => 0,
 				'MIN_%POS' => 0, 'MAX_%POS' =>101,
 				'MIN_BPGAP' =>0, 'MAX_BPGAP'=> 999999999999,
 				'MIN_%GAP' => 0, 'MAX_%GAP' =>101,
 				MIN_BPALIGN=> 0,       				###this matches plus mismatches 
 				MIN_FRACBPMATCH => 0, MAX_FRACBPMATCH =>1.01, # matches / (matches + mismatches)	
 				MIN_SUMBPALIGN=>0, MAX_SUMBPALIGN => 99999999999999,
 				NAME_TYPE => 'VERSION' , 
 		@_ ); 
		
 				
		die "Expecting a glob to the filehandle e.g. FH => \\*IN\n" if $args{'FILEHANDLE'} !~ /GLOB/;
		my $FH = $args{'FILEHANDLE'};  #too long type all the time
		my %d=();
		my @sbjct=();
		my %sbj=();


		my $line = '';	
		$line = <$FH> until $line=~ /^\#\:lav/ || eof $FH;
		if (eof $FH ) { 
			$d{'error'} = 'No lav title found';
            # DG Aug 2016
            print "No lav title found, returning\n";
            # end DG Aug 2016
			return \%d;
		} 
		#START PARSING SINGLE QUERY VS SINGLE SUBJECT WITH POSSIBLE REVERSE SUBJECT#
		print "$d{'errror'}\n";
		$line =<$FH> until $line=~ /d\s+\{/; #find d {
		my $lav_d='';
		$line =<$FH>;
		until ( $line=~/^\}/) {
			$lav_d.= $line;
			$line=<$FH>;
		}
		$d{'d'}=$lav_d;
		#print "$lav_d\n";
		#$d{'qname'}=$1;
		#if ($args{'NAME_TYPE'} eq 'VERSION') {
		#	$d{'qname'} = $1 if $d{'qname'} =~/\|([A-Z]+[_0-9.]+)\|/;
		#	$d{'qname'} =~ s/\.encode$//;
		#}

        # the following used to cause an infinite loop in the case
        # in which there was no s stanza in the file.  David Gordon,
        # Aug 2016
    
#		$line =<$FH> until $line=~ /s\s+\{/; #find h {

    

        until( $line=~ /s\s+\{/ ) {
            $line =<$FH>;
            if ( eof( $FH ) ) {
                print "no s { stanza in file so terminating\n";
                $d{'error'} = 'No s stanza found';

                return;                                                                    }
            print "looking for s { and found: " . $line
        }                                                                      

		####GET QUERY NAME AND LENGTH
		my $s1 =<$FH>;
		$s1=~/"(.*)"/;
		$d{'qname'}=$1;
		$d{'qname'}=$1 if $d{'qname'}=~/^.*\/(\S+)/ ;
		if ($args{'NAME_TYPE'} eq 'VERSION') {
			$d{'qname'} = $1 if $d{'qname'} =~/^([A-Z]+[_0-9]+\.[0-9_]+)/;
			$d{'qname'} =~ s/\.encode$//;
			$d{'qname'} =~ s/\.fugu$//;
		}

        # s stanza format changed from this with webb_self:
        # "/net/eichler/vol18/dgordon/wgac/hg19/fugu2/chr10_025.fugu" 1 183597
        # to this with lastz:
        # "fugu2/chr10_025.fugu" 1 183597 0 1
        # fixed David Gordon, Aug 25, 2016

		#print $d{'qname'},"\n";
		#$s1 =~ /(\d+)$/;


		$s1 =~ /\s+\d+\s+(\d+)\s+\d+\s+\d+$/ || die "s stanza didn't match pattern like  \"fugu2/chr10_025.fugu\" 1 183597 0 1\n";
		$d{'qlen'} =$1;
		#print "$d{'qlen'}\n";
		####GET SUBJECT NAME AND LENGTH
		my $s2 =<$FH>;
		$s2=~/"(.*)"/;
		$sbj{'name'}=$1;
		$sbj{'orient'}='F'; $sbj{'orient'}='R' if $sbj{'name'} =~/\-$/;
		print "SUB2:$sbj{'orient'}\n";
		$sbj{'name'}=$1 if $sbj{'name'}=~/^.*\/(\S+)/ ;
		if ($args{'NAME_TYPE'} eq 'VERSION') {
			$sbj{'name'} = $1 if $sbj{'name'} =~/^([A-Z]+[_0-9]+\.[0-9_]+)/;
			$sbj{'name'} =~ s/\.encode$//;
            # added to handle case in which there is only
            # 1 stanza and it is reverse complemented (David Gordon, Aug 26, 2016 )
			#$sbj{'name'} =~ s/\.fugu$//;
			$sbj{'name'} =~ s/\.fugu-?$//;
		}

        # s stanza format changed from this from webb_self:
        # "/net/eichler/vol18/dgordon/wgac/hg19/fugu2/chr10_025.fugu" 1 183597
        # to this from lastz:
        # "fugu2/chr10_025.fugu" 1 183597 0 1
        # fixed David Gordon, Aug 25, 2016

        #$s2 =~ /(\d+)$/;
        $s2 =~ /\s+\d+\s+(\d+)\s+\d+\s+\d+$/ || die "s stanza didn't match pattern like  \"fugu2/chr10_025.fugu\" 1 183597 0 1\n";
		$sbj{'len'} =$1;
		#print "SLEN:$sbj{'len'}\n";
		
		my $lav_s= $s1 . $s2;
		$d{'s'}=$lav_s;
		#print "S:$d{'s'}\n";
		
		$line =<$FH> until $line=~ /h\s+\{/; #find h {
		$line =<$FH>;
		my $lav_h='';
		until ( $line=~/^\}/) {
			$lav_h.= $line;
			$line=<$FH>;
		}
		$d{'h'}=$lav_h;
		#print "H:$d{h}\n";
		$lav_h=~ s/\r/\n/mg;
		$lav_h=~ s/\n/ /mg;
		$lav_h=~ s/ +/ /mg;
		$lav_h=~ s/^\|//;
		my $count=0;
		while ($lav_h =~ /">(.*?)"/mg){
			$d{'qdefn'}=$1 if $count==0;
			$sbj{'defn'}=$1 if $count==1;
			die "$count is too high " if $count > 1;
			$count++;
		}
		chomp $d{'qdefn'};
		chomp $sbj{'defn'};
		#print "QD:$d{'qdefn'}\nSD:$sbj{'defn'}\n";
		my @hsp;
		$line=<$FH> until $line =~/^[sa]/ || eof $FH;
		
S_or_A_STANZA:
		while ($line=~ /^[sa]/) {

		   if ($line=~/^a/) {
		    	my %h=();
		    	$line=<$FH>;
		    	( $h{'score'} ) = $line =~/(\d+)/;
		    	$line=<$FH>;
		    	( $h{'qb'},$h{'sb'} ) = $line =~ /(\d+)\s+(\d+)/;
		    	$h{'sb'}=$sbj{'len'}-$h{'sb'}+1 if $sbj{'orient'} eq 'R';
		    	$line=<$FH>;
		    	( $h{'qe'},$h{'se'} ) = $line =~ /(\d+)\s+(\d+)/;
		    	$h{'se'}=$sbj{'len'}-$h{'se'}+1 if $sbj{'orient'} eq 'R';
		    	$line=<$FH>;
		    	my ($bpident,$bpgap,$bpalign)=(0,0,0);
		    	my ($lqe,$lse)=(0,0);

L_STANZA:
		    	while ( $line=~/l (\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ ) {
		    		my ($qb,$sb,$qe,$se,$perc)=($1,$2,$3,$4,$5);
		    		my $leng= $qe-$qb+1;
		    		##this is the best guess since he is doing alot of rounding##
		    		#print "$line\n";
		    		if ($perc==0 || $leng<1) {
		  				$d{'error'}.= $line;
		  				print "BAD ALIGNMENT $line ($leng) at line $.\n";

                        # changed DG July 2019 to handle s 0 blocks
                        while( 1 ) {
                            $line=<$FH>;
                            if ( eof( $FH ) ) {
                                print "ignoring this last bad alignment\n";
                                return;####
                            }
                            elsif ( $line=~/l (\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ ) {
                                next L_STANZA;
                            }
                            elsif ( $line =~ /^[sa]/) {
                                next S_or_A_STANZA;
                            }
                        }

                        # end of DG change
 		  			}
		    		my $ident=int ($leng*($perc+0.499999)/100 ); #reverse engineering#
		    		my $fix=int($ident/$leng *100+ 0.5);  #forward check#
		    		#son of a bitch#
		    		#print "$leng ($perc) => $ident  fix ($fix)\n";
		    		$bpident+=$ident;
		    		$bpalign+=$leng;
		    		$line=<$FH>;
		    		if ($lqe>0) {
		    			###calculate gap statistics
		    			my $qgap=$qb-$lqe-1;
		    			my $sgap=$sb-$lse-1;
		    			if ($qgap>0 and $sgap>0) {
		    				print "PRINT DOUBLE GAP ERROR\n" ;
		    			}
		    			$bpgap+=$qgap + $sgap;
		    		}
		    		($lqe,$lse)=($qe,$se);
		    	}
		    	next if $bpalign <=0;
		    	my $bptot=$bpalign+$bpgap;
		    	my $fracbpalign=int($bpident/$bpalign*1000000)/1000000;
		    	my $percgap=int($bpgap/$bptot*10000)/100;
		    	my $percident=int($bpident/$bptot*10000)/100;
		    	( $h{'bpident'},$h{'sizealign'}, $h{'%ident'}, $h{'bpgap'})=($bpident,$bptot,$percident,$bpgap);
		    	( $h{'%gap'}, $h{'bpalign'}, $h{'fracbpmatch'}) = ($percgap,$bpalign,$fracbpalign);
		    	#print "OVERALL:$h{'score'} $h{qb}($h{sb}) to $h{qe}($h{se}) A: $h{'bpident'}/$h{'bpalign'} ($h{'fracbpmatch'})"
		    	#	." I:$h{bpident}/$h{'sizealign'} ($h{'%ident'})"
		    	#	." G:$h{bpgap}/$h{'sizealign'} ($h{'%gap'})\n";
		    	
				############################################ 
				#cue to next record
				$line=<$FH> until $line =~/^[sa]/ || eof $FH; 
				###SKIP AN INDVIDUAL HSP WITH A NEXT ####
				
				next if $h{'bpident'} < $args{'MIN_BPIDENT'};
				next if $h{'bpalign'} < $args{'MIN_BPALIGN'};
				
				next if $h{'fracbpmatch'} < $args{'MIN_FRACBPMATCH'};
				next if $h{'fracbpmatch'} > $args{'MAX_FRACBPMATCH'};
				#print "DEciding ...\n";
				next if $h{'sizealign'} < $args{'MIN_SIZEALIGN'};
				next if $h{'%gap'} > $args{'MAX_%GAP'};
				next if $h{'%gap'} < $args{'MIN_%GAP'};
				######what to do if skipsame is false############
				if ($args{'SKIP_SAME'} == $false) {
				   #remove identical hits
				next if ( $d{'qname'} eq $sbj{'name'} ) && ($h{'qb'} == $h{'sb'}) && $args{'SKIP_IDENT_HSP'};
				#remove mirrors in for self hits#
				my $skip= $false;
				foreach my $x (@hsp) {
					$skip=$true if  $h{'qb'}==$$x{'sb'} && $h{'qe'}==$$x{'se'} && $h{'sb'}==$$x{'qb'} && $h{'se'}==$$x{'qe'};
					$skip=$true if  $h{'qb'}==$$x{'se'} && $h{'qe'}==$$x{'sb'} && $h{'sb'}==$$x{'qe'} && $h{'se'}==$$x{'qb'};
					}
				next if $skip;
				}
				#print "KEEPER\n";
				push @hsp, \%h;		    	
		    	next;
		    } elsif ($line=~/^s/) {
		    		#print "Reverse Complement...";
		    		$line=<$FH>; $line=<$FH>;
		    		if ($line=~ /\-"/) {
		    			$sbj{'orient'}='R';
		    			print "ORIENT$sbj{'orient'} $line\n";
		    			$line=<$FH> until $line =~/^[sa]/ || eof $FH;
		    		} else {
		    			die "I don't know why there are two files together\n";
		    		}
		    } else {
		    	die "Parsing error\n";
		    }
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

		###PROCESSING BASED ON THE SUBJECT AND HSP TOTALS ##########	
		$sbj{'sumbpalign'}=0;
		foreach my $h (@hsp) { $sbj{'sumbpalign'}+=$$h{'bpalign'}; }
		next if $sbj{'sumbpalign'} < $args{'MIN_SUMBPALIGN'};
		next if $sbj{'smubpalign'} > $args{'MAX_SUMBPALIGN'};

		next if  ( $d{'qname'} eq $sbj{'name'} ) && $args{'SKIP_SELF'};
		next if  $args{'SKIP_SUBJECTS'} && $args{'SKIP_SUBJECTS'}{$sbj{'name'}};

		push @sbjct, \%sbj;

		#####PROCESSING BASED ON HAVING ALL SUBJECTS AVAILABLE ####
		##########################################################
		@{ $d{'sbjct'} }=@sbjct;
		#$line =<$FH> until $line=~/^BLAST/ || eof $FH;
		return \%d;



}
	


sub parse_query_format_pairwise {
	my %args=(FILEHANDLE=>\*STDOUT, QUERY=>'', @_);
	my $fh= $args{'FILEHANDLE'};
	foreach my $s ( @{ $args{'QUERY'}{'sbjct'} }) {
		my $q_same_s=$false;
		print "$$s{'name'} ";
		if ( $args{'QUERY'}{'qname'} eq $$s{'name'}  ) {
			foreach my $h  ( @{$$s{'hsps'}}  ) {
				#print  "\n$args{'QUERY'}{qname}\t$$h{qb}--$$h{qe}\t$args{'QUERY'}{'qlen'}";
				#print "\t$$s{name}\t$$h{sb}--$$h{se}\t$$s{len}\n";
				
				if ($$h{'sb'} < $$h{'qb'} || $$h{'se'} < $$h{'qb'} ) {
					if ( $$h{'sb'} < $$h{'se'} ) {
						##subject forward##
						($$h{'qb'},$$h{'qe'},$$h{'sb'},$$h{'se'})=($$h{'sb'},$$h{'se'},$$h{'qb'},$$h{'qe'});
					} else {
						($$h{'qe'},$$h{'qb'},$$h{'se'},$$h{'sb'})=($$h{'sb'},$$h{'se'},$$h{'qb'},$$h{'qe'});
					}
				}
				#print  "$args{'QUERY'}{qname}\t$$h{qb}--$$h{qe}\t$args{'QUERY'}{'qlen'}";
				#print "\t$$s{name}\t$$h{sb}--$$h{se}\t$$s{len}\n";
				#my $pause=<STDIN>;

			}
		}
		@{ $$s{'hsps'}   } = sort { $$a{'qb'} <=> $$b{'qb'} } (@{ $$s{'hsps'}   } );
		
		foreach my $h ( @{ $$s{'hsps'}   }) {
			print $fh "$args{'QUERY'}{qname}\t$$h{qb}\t$$h{qe}\t$args{'QUERY'}{'qlen'}";
			print $fh "\t$$s{name}\t$$h{sb}\t$$h{se}\t$$s{len}";
			print $fh "\t$$h{'fracbpmatch'}\t$$h{'bpalign'}\t$$h{'sizealign'}";
			print $fh "\t$$h{'score'}";
			print $fh "\t",$args{'QUERY'}{'qdefn'};
			print $fh "\t","$$s{'defn'}";
			print $fh "\n";
		}
	}
	print "\n";
}

	
	
