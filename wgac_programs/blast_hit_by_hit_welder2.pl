#!/usr/bin/perl


###VERSION 0.000807 alpha####

##000807 added htgs draft grouping and better output and generally finished program#####
##CREATED ON 24-07-00

use strict 'vars';
use Getopt::Std;
use Data::Dumper;

use vars qw($program $version);
use vars qw ($true $false);
use vars qw ($smallergap $largergap $overlap);
use vars qw ($opt_i $opt_s $opt_l $opt_g $opt_o $opt_s);
use vars qw (%h);

($true,$false)=(1,0);

$program = "GB4_clean_current_versions.pl";
$version="2.001001";

if (!defined $ARGV[0]) {
print "USAGE blast_hit_by_hit_welder.pl -i table1:table2 -s [min gap size] -l [max gap size]
REQUIRED ARGUMENTS
-i [table paths]  colon delimited paths of tables to be combined before welding them together
            **these tables will be sorted and require for self hits that seq1 is lower number**
            **also requires that between any given 2 sequences s1 and s2 remain the same for all hits between them**
-g [smallergap:largergap:overlap] small:large:overlap
#-c [maximum finished gap] to htgs for htgs clustering
-l [integer]   the maximum size for the larger gap (default 5000 bp)
-s [interger] start number 
-o [path] base output path
";
exit;
}


getopts('i:s:l:g:o:');
$opt_i || die "Please use -i with : delimits to  designate blast hit tables to weld!\n";

my @tables = split ":", $opt_i;
print "TABLES:",@tables,"\n";


die "Use : colons to delimit -g  (e.g. 100:5000:25)\n" if $opt_g && $opt_g !~ /:/;
($smallergap,$largergap,$overlap)=split ":", $opt_g;
$smallergap = 10 if $smallergap eq '';
$largergap = 2000 if $largergap eq '';
$overlap = 10 if $overlap eq '';

print "GAP PARAMETERS ($smallergap)::($largergap)::($overlap)\n";

for (my $i=0;$i<@tables;$i++) {
	&blast_table_sorter(-table=> $tables[$i]);
	$tables[$i].='.sort';
};

&blast_table_combiner(-tables => \@tables, -out => "$opt_o.combocull");
&blast_welder (-table=>"$opt_o.combocull", -small=>$smallergap, -large=>$largergap, 
			-out=> "$opt_o",-overlap=>$overlap);




sub blast_table_sorter {
	##This subroutine will sort tabdelimited table#
	my $sub_name='table_sorter';
	##set defaults and check for parameter errors
	my @defaults= (-table=>'', -out => '');
	my %defaults=@defaults;
	my %args= (@defaults,@_);
	foreach my $k (keys %args) {
		if (!defined $defaults{$k} ) {
			print "\n$k is an invalid parameter for subroutine $sub_name\n";
			print "  VALID PARAMETERS for $sub_name are: ";
			print join ("':'",keys %defaults),"\n";
			die "Execution halted\n";
		}
	}
	my @rows=();
	my $t=$args{'-table'};
	my $out=$args{'-out'};
	open (IN, $t) || die "Can't open input table to sort ($t)!\n";
	my $header=<IN>;
	#print $header,"\n";
	while (<IN>) {
		s/\r\n//;
		chomp;
		my @c = split "\t";
		push @rows, \@c;
	}
	close IN;
	@rows=sort {$$a[0] cmp $$b[0] or $$a[4] cmp $$b[4] or $$a[1] <=> $$b[1]} @rows;
	
	$out ||= "$t.sort";
	open (OUT,">$out") || die "Can't open sorted output table ($out)!\n";
	print OUT $header;
	foreach (@rows) {
		print OUT join ("\t", @$_), "\n";
	}
	close OUT;
}	
	
	
	
sub blast_table_combiner {
	##This subroutine will combine sorted tables#
	my $sub_name='table_combiner';
	
	##set defaults and check for parameter errors
	my @defaults= (-tables=>'', -out => '', -min_bpalign=> 100, -min_fracbpalign=>0.90);
	my %defaults=@defaults;
	my %args= (@defaults,@_);
	foreach my $k (keys %args) {
		if (!defined $defaults{$k} ) {
			print "\n$k is an invalid parameter for subroutine $sub_name\n";
			print "  VALID PARAMETERS for $sub_name are: ";
			print join ("':'",keys %defaults),"\n";
			die "Execution halted\n";
			
		}
	}
	
	my @tables = @{$args{'-tables'}};
	###MAIN ALGORITHM###
	my %line;
	my $header='';
	foreach my $t (@tables) {
		open ($t,"$t") || die "Can't open ($t) for reading\n";
		$header =<$t>;      #grab the header#
		$line{$t}=<$t>;  #prime the array#
	}
	open (OUT, ">$args{'-out'}" ) || die "Can't open final output table ($args{-out})!\n";
	print OUT $header;
	
	my $all_eof=$false;
	while ($all_eof==$false) {
		my ($first,$seqf,$beginf);
		foreach my $t (@tables) {
			next if $line{$t} eq "EOF99";
			my ($seq, $begin)=split "\t",$line{$t};
			if (! defined $first) {
				$first=$t; $seqf=$seq; $beginf=$begin;
			}
			if ($seq lt $seqf || ($seq eq $seqf && $begin < $beginf) ) {
				$first=$t; $seqf=$seq; $beginf=$begin;
			}
			
		}
		#print "$seqf $beginf ($first)\n";
		my @c=split "\t",$line{$first};
		chomp @c;
		if ( $c[27] >=$args{-min_fracbpalign} && $c[22]>=$args{-bpalign}) {
			print OUT $line{$first} ;
		}
		if (!eof $first) {
				$line{$first} = <$first>;
		} else {
			$line{$first} = "EOF99";
			$all_eof=$true;
			foreach my $t2 (@tables) {
				if ($line{$t2} ne "EOF99") {
					$all_eof=$false;
					last;
				}
			}
		}	
	}
	close OUT;
}



sub blast_welder {
	print "WELDING\n";
	##This subroutine will combine sorted lists##
	my $sub_name='blast_welder';
	##set defaults and check for parameter errors
	my @defaults= (-table=>'', -small=> 1, -large => 1000, -overlap=> 1, -out => '', -dpattern=>'P\:[01]\;',
						 -dmax=>1000, -doverlap=>1, -min_bpalign=>500);
	my %defaults=@defaults;
	my %args= (@defaults,@_);
	foreach my $k (keys %args) {
		if (!defined $defaults{$k} ) {
			print "\n$k is an invalid parameter for subroutine $sub_name\n";
			print "  VALID PARAMETERS for $sub_name are: ";
			print join ("':'",keys %defaults),"\n";
			die "Execution halted\n";
			
		}
	}
	print "SMALL:$args{-small} LARGE:$args{-large} OVERLAP $args{-overlap}\n";
	
	open (TABLE,$args{'-table'}) || die "Can't open table ($args{'-table'})\n";
	my @t=<TABLE>;
	close TABLE;
	
	###MAIN ALGORITHM TO DETERMINE THE FINISHED SEQUENCE WELDS###
	my $weld=0;
	my $weldcount=0;
	my @w;
	print "FINISHED WELDING ($args{-small}:$args{-large}:$args{-overlap})\n";
	for (my $i=1;$i<@t;$i++) {
		my @ai=split "\t", $t[$i];
		my $orienti=$ai[6]-$ai[5]; if ($orienti<0) {$orienti='-';} else { $orienti='+';}
		my $phase1="F";   $phase1="D" if $ai[12] =~ /$args{'-dpattern'}/;
		my $phase2 = "F"; $phase2="D" if $ai[13] =~ /$args{'-dpattern'}/;
		next if $phase2 eq 'D';

		for ( my$j=$i+1; $j<@t; $j++) {
			my @aj=split "\t",$t[$j];
			my $orientj= $aj[6]-$aj[5];
			if ($orientj<0) {$orientj='-';} else { $orientj='+';}
			###lasts###
			last if $aj[0] ne $ai[0]; #seq names at col 0;
			my $seq1diff=$aj[1]-$ai[2]-1;
			#print "($i:$ai[1])S1: $ai[2] -- $aj[1]  => $seq1diff\n";
			if (  $seq1diff > $args{'-large'} ) {
				#print "TOO LARGE\n";
				next;
			}
			if ( $seq1diff < -$args{'-overlap'}) {
				#print "TOO MUCH OVERLAP\n";
				next;
			}
			if ( $orienti ne $orientj ) { 
				#print "BAD ORIENTATION\ $orienti ne $orientj \n ";
				next;
			}
			if ( $aj[4] ne $aj[4] ) {next;}  #seq names at col 4 #keep looking for seq2 name#;
			my $seq2diff=$aj[5]-$ai[6]-1;
			if ($orienti eq '-') {
				$seq2diff=$ai[6]-$aj[5]-1 ;	
			}
			#print "  $orienti:$orientj  $ai[6] -- $aj[5]  ==> $seq2diff\n";
			
			if ($seq2diff < -$args{'-overlap'} ) {
				#print " SEQ2 TOO MUCH OVERLAP\n";
				next;
			}
			if ($seq2diff > $args{'-large'}) {
				#print "SEQ2 TOO LARGE\n";
				next;
			}
			if ($seq1diff > $args{'-small'} && $seq2diff > $args{'-small'} ) {
				#print "BOTH GAPS TOO BIG\n";
				next;
			}
			if ($orienti eq '+' && $ai[5] < $aj[2] && ($ai[0] eq $ai[4] )) { 
				print "WARNING this next one self overlaps *****\n";
				next;
			}
			if ($orienti eq '-' && $aj[6] < $aj[2] && ($ai[0] eq $ai[4])) {  ###this should be impossible but added for completeness###
				print "WARNING this next one self overlaps  *****\n";
				next;
			}			
			if (defined $w[$j]) {
				print "WARNING PROMISCUS J PIECE ALREADY DEFINED ($w[$j])\n";
			}
			if (!defined $w[$i] ) {
				$w[$i]= "WELD$phase1$phase2" . ++$weld;
			} 
			$w[$j]=$w[$i];
			
			print " $w[$i]---\n";
			print "  $phase1$i$orienti:$phase2$j$orientj  12 $ai[0] $ai[1] $ai[2] ($aj[0] $aj[1] $aj[2])  [$seq1diff]\n";
			print  "    56  $ai[4] $ai[5] $ai[6] ($aj[4] $aj[5] $aj[6])  [$seq2diff]\n";
			$weldcount++;
			#my $pause=<STDIN>;
		}
	}
	
	my $group=10000;
	my @h=();
	print "DRAFT GROUPING...  ($args{'-dpattern'})($args{-dmax})($args{-doverlap})\n";
	ILOOP:for (my $i=1;$i<@t;$i++) {
		my @ai=split "\t", $t[$i];
		my $phase1="F"; $phase1="D" if $ai[12] =~ /$args{'-dpattern'}/;
		my $phase2="F"; $phase2="D" if $ai[13] =~ /$args{'-dpattern'}/;
		next if $phase1 eq 'D';
		next if $phase2 eq 'F';
		#print "PHASE $phase1 $phase2\n";
		JLOOP:for ( my$j=$i+1; $j<@t; $j++) {
			my @aj=split "\t",$t[$j];
			my $seq1diff = $aj[1]-$ai[2]-1;
			next ILOOP if  $aj[4] ne $ai[4];
			next JLOOP if  $seq1diff > $args{'-dmax'};
			next JLOOP if  $seq1diff < -$args{'-doverlap'};
			next JLOOP if defined $h[$i] && defined $h[$j] && $h[$i] eq $h[$j];
			
			if (defined $h[$j]) {
				if (!defined $h[$i]) {
					$h[$i]=$h[$j];
					next JLOOP;
				}
				print "WARNING PROMISCUS I$i J$j PIECE ALREADY DEFINED ($h[$j])\n";
				
				print "$h[$i]--$h[$j]-\n";
				print "$phase1$i:$phase2$j  12 $ai[0] $ai[1] $ai[2] ($aj[0] $aj[1] $aj[2]) ($ai[4]:$aj[4])  [$seq1diff]\n";
				next JLOOP;
			}
			if (!defined $h[$i] ) {
				$h[$i]= "GROUPFD" . ++$group;
			} 
			$h[$j]=$h[$i];
			
			#print "$h[$i]--$h[$j]-\n";
			#print "$phase1$i:$phase2$j  12 $ai[0] $ai[1] $ai[2] ($aj[0] $aj[1] $aj[2]) ($ai[4]:$aj[4])  [$seq1diff]\n";

			
			
		}
	}
		
	###Collapse w and H to W###	
	ILOOP:for (my $i=1;$i<@t;$i++) {	$w[$i]=$h[$i] if defined $h[$i]; }
 	#save t for outputing all alignments
	my @t_save=@t;
	my @w_save=@w;
	##########################################################################################
	##########MEGE ALIGNMENT FILES#########
	print "MERGING ALIGNMENT FILES...\n";
	my %welds;
	ILOOP:for (my $i=1;$i<@t;$i++) {
		#die if $i==22;
		next if  ! defined $w[$i];
		my @ai=split "\t", $t[$i];
		my $orienti=$ai[6]-$ai[5]; if ($orienti<0) {$orienti='-';} else { $orienti='+';}
		print "  $w[$i]---\n";
		print "   $i: $ai[0] $ai[1] $ai[2]    $ai[4] $ai[5] $ai[6]  ($w[$i])\n";
		JLOOP:for (my $j=$i+1;$j<@t;$j++) {
			my @aj=split "\t", $t[$j];
			my $gap1= $aj[1]-$ai[2]-1;
			my $gap2=$aj[5]-$ai[6]-1;   $gap2=$ai[6]-$aj[5]-1 if $orienti eq '-';
			next ILOOP if $ai[0] ne $aj[0];
			next ILOOP if $ai[4] ne $aj[4];
			next ILOOP if $gap1 > 20000;
			next JLOOP if $w[$j] ne $w[$i];
			print "   *$gap1*$j: $aj[0] $aj[1] $aj[2]   $aj[4] $aj[5] $aj[6]  ($w[$j])\n";
			####combine the two
			#c[0] nothing
			#c[1,2] 
			$ai[2] = $aj[2];  #the longer of the two
			#c[3,4] nothing
			#c[5,6]
			if ($w[$i]=~ /WELD/ ) {
				if ($orienti eq '+') { $ai[6]=$aj[6]; } else { $ai[6]=$aj[6]; }
			} else {
				my ($l,$h) = ($aj[5],$aj[6]);
				($l,$h)=($h,$l) if $l>$h;
				if ($orienti eq '+') {
					$ai[5]=$l if $l<$ai[5];
					$ai[6]=$h if $h>$ai[6];
				} else {
					$ai[5]=$h if $h>$ai[5];
					$ai[6]=$l if $l<$ai[6];
				}
			}
			#c[7] nothing
			#c[8-15]  blast crap ignore
			for (my $k=8;$k<=11; $k++) { $ai[$k]='NA' ; }
			#c[16,17] append ...
			$ai[16].=":$aj[16]";
			$aj[17].=":$aj[17]";
			#approximate gap caculations#
			my $gapspace;
			if ($gap1<$gap2) { $gapspace=$gap2-$gap1; } else { $gapspace=$gap1-$gap2;}
			#c18 nothgin
			for (my $k=19;$k<=26; $k++) { $ai[$k]+=$aj[$k] ; }
			
			$ai[20]+=1;
			$ai[27]= $ai[23]/($ai[22]);
			$ai[28]='NA';
			$ai[29]= $ai[23]/($ai[22]+$ai[20]);  #sim with indels
			$ai[30]='NA';
			($ai[31],$ai[32]) = &k_jukes_cantor($ai[27],$ai[22]);
			($ai[33],$ai[34]) = &k_kimura($ai[26],$ai[25],$ai[22]); #transitions,transversion,base_spaces
			
			for (my $k=35;$k<=37; $k++) { $ai[$k]='NA' ; }
			for (my $k=38;$k<=62; $k++) { $ai[$k]+=$aj[$k] ; }
			$t[$i]= join( "\t", @ai);
			$welds{$w[$i]}{'bpalign'}=$ai[22];
			splice (@t,$j,1); 
			splice (@w,$j,1);
			$i--;
			next ILOOP;
			
			
		}
		
	}
	
	




	####WRITE OUT FILE WITH JOINED WELDS NUMBERS IN LAST COLUMN#####
	print "OUTPUTING JOINS...\n";
	open (OUT, ">$args{-out}.join.all") || die "Can't open out table ($args{-out}.joinall) for writing!\n";
	open (OUTW, ">$args{-out}.join.weld") || die "Can't open out table ($args{-out}.justweld) \n";
	open (OUTG, ">$args{-out}.join.group") || die "Can't open out table ($args{-out}.justgroup)\n";
	open (OUTNJ, ">$args{-out}.join.nojoin") || die "Can't open out table ($args{-ou8t}.nojoin)\n";
	
	chomp @t;
	print OUT "$t[0]\tWELD\n";
	print OUTW "$t[0]\tWELD\n";
	print OUTG "$t[0]\tWELD\n";
	print OUTNJ "$t[0]\tWELD\n";
	for (my $i=1;$i<@t;$i++) {
		my @c=split "\t",$t[$i];
		chomp @c;
		next if  ($c[22] < $args{'-min_bpalign'});

		print OUT $t[$i],"\t";
		if (defined $w[$i]) {
			print OUT "$w[$i]\t\n"; 
		} else { 
			print OUT "\n";
		}
		if (defined $w[$i] ) {
			if ($w[$i]=~/WELD/) {
				print OUTW $t[$i],"\t$w[$i]\n";
			} else {
				print OUTG $t[$i],"\t$w[$i]\n"
			}
		} else {
				print OUTNJ $t[$i],"\t\n";
		}
			
	}
	close OUT;
	close OUTW;
	close OUTG;
	close OUTNJ;
	
	@t=@t_save; 
	@w=@w_save;
	####WRITE OUT FILE WITH TRAILING WELD NUMBERS IN LAST COLUMN#####
	print "OUTPUTING PIECES...\n";
	open (OUTNW, ">$args{-out}.pieces.weld") || die "Can't open out table ($args{-out}.pieces.weld)!\n";
	open (OUT, ">$args{-out}.pieces.all") || die "Can't open out table ($args{-out}.pieces.all) for writing!\n";
	chomp @t;
	print OUT "$t[0]WELD\n";
	print OUTNW "$t[0]WELD\n";
	for (my $i=1;$i<@t;$i++) {
		my @c=split "\t",$t[$i];
		chomp @c;
		next if ( $c[22]<$args{'-min_bpalign'} && defined $welds{$w[$i]} && $welds{$w[$i]}{'bpalign'}<$args{'-min_bpalign'});
		if ($w[$i] !~ /WELD/ ) {
			print OUTNW $t[$i],"\t";
			if (defined $w[$i]) {print OUTNW "$w[$i]\t$welds{$w[$i]}{'bpalign'}\n"; } else { print OUTNW "NO\t$c[22]\n";}
		}
		print OUT $t[$i],"\t";
		if (defined $w[$i]) {print OUT "$w[$i]\t$welds{$w[$i]}{'bpalign'}\n"; } else { print OUT "NO\t$c[22]\n";}

	}
	close OUT;	
	print "WELDS:$weld JOINTS:$weldcount\n";

}

sub k_jukes_cantor {
	my $frac_sim=shift; #frac similarity mismatches/base spaces(aligned bases no indel spaces)
	my $p=1-$frac_sim;
	my $base_spaces=shift; #basespaces
	#print "$p $base_spaces\n";
	my $k_jukes_cantor= -0.75*log(1-4/3*($p));
	my $SE_jukes_cantor='NA'; 
	$SE_jukes_cantor=((1-$p)*$p/($base_spaces*(1-4*(1-$p)/3)**2))**0.5 if $base_spaces;
	return $k_jukes_cantor, $SE_jukes_cantor;
}

sub k_kimura { 
	my ($transitions,$transversions,$base_spaces)=@_;
	my $p =$transitions/$base_spaces;
	my $q = $transversions/$base_spaces;
	my $a=1/(1-2*$p-$q);
	my $b=1/(1-2*$q);
	my $k_kimura= 0.5 * (log $a) + 0.25*( log $b);
	my $SE_kimura= ( ($a**2*$p+ (($a+$b)/2)**2*$q - ($a*$p + ($a+$b)/2*$q )**2)/$base_spaces )**0.5;
	return $k_kimura, $SE_kimura;
}