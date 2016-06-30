#!/usr/bin/perl



###welding together post alignments####

use strict 'vars';
use Getopt::Std;
use Data::Dumper;

use vars qw($program $version);
use vars qw ($true $false);
use vars qw ($smallergap $largergap $overlap);
use vars qw ($opt_i $opt_s $opt_l $opt_g $opt_o $opt_s);
use vars qw (%h);

($true,$false)=(1,0);

$program = "";
$version="3.020926";

#02-09-26 fixed the problem with out being screwed up for row #s  *.fixed and *.fixed.trim are identical
###I really have no idea what the parameters are for the merging

if (!defined $ARGV[0]) {
print "
REQUIRED ARGUMENTS
-i [table paths]  colon delimited paths of tables to be combined before welding them together
            **these tables will be sorted and require for self hits that seq1 is lower number**
            **also requires that between any given 2 sequences s1 and s2 remain the same for all hits between them**
			This program is for merging data
			very poorly annotated.
";
exit;
}


getopts('i:s:l:g:o:f:');
$opt_i || die "Please use -i with : delimits to  designate blast hit tables to weld!\n";
$opt_o || ($opt_o=$opt_i);

my @tables = split ":", $opt_i;
print "TABLES:",@tables,"\n";


die "Use : colons to delimit -g  (e.g. 100:5000:25)\n" if $opt_g && $opt_g !~ /:/;
($smallergap,$largergap,$overlap)=split ":", $opt_g;
$smallergap = 0 if $smallergap eq '';
$largergap = 0 if $largergap eq '';
$overlap = 10 if $overlap eq '';

print "GAP PARAMETERS ($smallergap)::($largergap)::($overlap)\n";

print "SORTING TABLES\n";
for (my $i=0;$i<@tables;$i++) {
	&blast_table_sorter(-table=> $tables[$i]);
	$tables[$i].='.sort';
};
print "COMBINING TABLES\n";
&blast_table_combiner(-tables => \@tables, -out => "$opt_o.combocull");
print "RUNNING BLAST FIXER\n";
&blast_fixer (-table=>"$opt_o.combocull", -small=>$smallergap, -large=>$largergap, 
			-out=> "$opt_o.fixed",-overlap=>$overlap);

&trim_format ( -in=> "$opt_o.fixed" , -out => "$opt_o.fixed.trim"  );

sub trim_format {
	my $sub_name='trim_format';
	##set defaults and check for parameter errors
	my @defaults= (-in=>'', -out => '');
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
	open (IN, "$args{-in}") || die "Can't read ($args{-in})!\n";

	open (OUT, ">$args{'-out'}") || die "Can't write ($args{'-out'})!\n";
	my $head=<IN>;
	$head =~ s/FILE\tBEGIN\t.*//;
	print OUT "$head";
	my $count=0;
	while (<IN>) {
		#print;
		$count++;
		my $row = substr ("0000000$count",-7);
		s/\trow\d+.*/\trow$row/;
		print  OUT ;
		print "$count..\n" if $count % 5000 ==0;
		#my $pause=<STDIN>;
	}
	close IN;
	close OUT;	
}


sub blast_table_sorter {
	##This subroutine will sort tabdelimited table##
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
	@rows=sort {$$a[0] cmp $$b[0] or $$a[4] cmp $$b[4] or $$a[1] <=> $$b[1] or $$b[2]<=> $$a[2] } @rows;
	
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
	my @defaults= (-tables=>'', -out => '', -min_bpalign=> 100, -min_fracbpalign=>0.88);
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
		if ( $c[8] >=$args{-min_fracbpalign} && ($c[2]-$c[1])>=$args{-bpalign}) {
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



sub blast_fixer {
	print "WELDING\n";
	##THIS SUBROUTINE OVERLAPPING PAIRWISE AFTER TRIMMING #####
	###########################################################
	##IT DOES NOT USE DMAX OR DOVERLAP
	###########################################################
	my $sub_name='blast_fixer';
	##set defaults and check for parameter errors
	my @defaults= (-table=>'', -small=> 1, -large => 1000, -overlap=> 1, -out => "$opt_o.fix", -dpattern=>'P\:[01]\;',
						 -dmax=>1000, -doverlap=>1, -min_bpalign=>1000);
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
	print "SMALL:$args{-small}\nLARGE:$args{-large}\nOVERLAP $args{-overlap}\n";
	print "DMAX:$args{-dmax}\nDOVERLAP:$args{-doverlap}\nMINBPALIGN $args{-min_bpalign}\n";
	
	open (TABLE,$args{'-table'}) || die "Can't open table ($args{'-table'})\n";
	my @t=<TABLE>;
	close TABLE;
	
	open (OUT, ">$args{'-out'}") || die "Can't write to ($args{'-out'})!\n";
	
	###MAIN ALGORITHM TO DETERMINE THE FINISHED SEQUENCE WELDS###
	my $weld=0;
	my $weldcount=0;
	my @w;
	ILOOP:for (my $i=1;$i<@t;$i++) {
		my @ai=split "\t", $t[$i];
		my $orienti=$ai[6]-$ai[5]; if ($orienti<0) {$orienti='-';} else { $orienti='+';}

		for ( my$j=$i+1; $j<@t; $j++) {
			my @aj=split "\t",$t[$j];

			###Simple checks#####
			last if $aj[0] ne $ai[0]; #seq names at col 0;
			last if $ai[4] ne $aj[4];
			next if $ai[0] ne $aj[0]; 
			next if $aj[1]- $ai[2] > 1000;

			my $orientj= $aj[6]-$aj[5];
			my $si12= $ai[2]-$ai[1]+1;
			my $sj12= $aj[2]-$aj[1]+1;
			my $si56= abs($ai[5]-$ai[6]);
			my $sj56 =abs($aj[5]-$aj[6]);

			if ($orientj<0) {$orientj='-';} else { $orientj='+';}
			###lasts###
			
			#next if $ai[0] eq $ai[4];  #if same skip why????


			#####skip if the join would be between different orienations####			
			if ( $orienti ne $orientj ) { 
				#print "BAD ORIENTATION\ $orienti ne $orientj \n ";
				next;
			}
			print "######$i)$orienti:$j) $ai[0] and $ai[4]\n";
			my $seq1diff;	
			my $contained1=0;
			if ( $ai[1]<= $aj[1]) {
				if ($ai[2]>=$aj[2]) {
					print "$ai[1]-----------------------------$ai[2]\n";
					print "         $aj[1]------$aj[2]\n";
					$seq1diff=$aj[1]-$ai[2];
					$contained1=1;
				}else {
				   print "$ai[1]----$ai[2] ..... $aj[1]----$aj[2]\n";
				   $seq1diff=$aj[1]-$ai[2];
				   
				}
			} else {
				die " $aj[1] is greater than $ai[1]--not sorted!\n";
			}
			print "##SEQ1DIFF $seq1diff ($contained1)\n";
			
			#plus orient#
			my $seq2diff=0;
			my $contained2=0;
			if ($orienti eq '+') {
				if ($ai[5]<=$aj[5]) {
					if ($ai[6]>=$aj[6]) {
						print "$ai[5]---------->----------->--------$ai[6]\n";
						print "         $aj[5]---->--$aj[6]\n";
						$seq2diff=$aj[5]-$ai[6];
						$contained2=1;
					}else {
						print "$ai[5]--->-$ai[6] ..... $aj[5]-->--$aj[6]\n";
						$seq2diff=$aj[5]-$ai[6];
					}
				} else {
					$seq2diff=99999999;
					print "$aj[5] is less then $ai[5] switchero\n";
					print "ai12 --> aj12\naj56--> ai56\n";
				}
			} else {
				#negative orientation#
				if ($ai[5] >= $aj[5]) {
					if ($ai[6]<=$aj[6]) {
						print "$ai[5]<--<----------<------<--------<------$ai[6]\n";
						print "          $aj[5]<-----<---$aj[6]\n";
						$seq2diff=$ai[6]-$aj[5];
						$contained2=1;
					} else {
						print "$ai[5]<--<---$ai[6] ... $aj[5]<--<---$aj[6] \n";
						$seq2diff=$ai[6]-$aj[5];
						
					}
				} else {
					print "$aj[5] is greater than $ai[5] (switcheroed)\n";
					print "$ai[5]<--<---$ai[6] wrong side $aj[5]<--<---$aj[6] \n";
				}
			}
			print "##SEQ2DIFF $seq2diff ($contained2)\n";
			#my $pause=<STDIN>;


			#print "12 $ai[0] $ai[1] $ai[2] ($aj[0] $aj[1] $aj[2])  [$seq1diff]\n";
			#print "   $si12      $sj12\n";
			#print "   $si56      $sj56\n";
			#print  "56  $ai[4] $ai[5] $ai[6] ($aj[4] $aj[5] $aj[6])  [$seq2diff]\n";
			next if  ($seq1diff >=-10 || $seq2diff>=-10)  && !($contained1 && $contained2) ;

			#my $pause=<STDIN>;



			my $min= $seq1diff;
			my $minspan=$si12;
			$minspan=$si56 if $si56 < $minspan;
			$minspan=$sj12 if $sj12 < $minspan;
			$minspan=$sj56 if $sj56 < $minspan;
			$min=$seq2diff if $seq2diff < $seq1diff;
			my $lognum=$minspan-abs($seq1diff-$seq2diff);
			$lognum =1 if $lognum < 1;
			my $max_gap = sqrt($lognum)  * 30;
			print "MINSPAN:$minspan ", abs($seq1diff-$seq2diff)," $min $max_gap\n";
			if ( $ai[1] eq $aj[1] && $ai[5] eq $aj[5]) {
						print "***BEGIN is EXACT MATCH\n";
			} elsif ($ai[2] eq $aj[2] && $ai[6] eq $aj[6]) {
						print "***END MATCH is EXACT MATCH\n";
			} elsif (abs($seq1diff -$seq2diff) < $max_gap ) {
				print "GOOD JOIN\n"  ;
				if ( $ai[0] eq $aj[0] ) {
					#extra careful melding overlaps if melding overlaps#
						print "NOT EXACTING skiiping\n";
						#my $pause=<STDIN>;
						next;
				}
			} else {
				print "TOO DIFFERENT\n";
				#my $pause=<STDIN>;
				next;
			}
			my ($nqb,$nqe,$nsb,$nse);
			if ($ai[1]<$aj[1]) { $nqb=$ai[1];} else {$nqb=$aj[1];}
			if ($ai[2]>$aj[2]) { $nqe=$ai[2];} else {$nqe=$aj[2];}
			if ($orienti eq '+' ) {
				if ($ai[5]<$aj[5]) { $nsb=$ai[5];} else {$nsb=$aj[5];}
				if ($ai[6]>$aj[6]) { $nse=$ai[6];} else {$nse=$aj[6];}
			} else {
				if ($ai[5]>$aj[5]) { $nsb=$ai[5];} else {$nsb=$aj[5]}
				if ($ai[6]<$aj[6]) { $nse=$ai[6];} else {$nse=$aj[6];}
			
			}
			print "NEW: ($nqb--$nqe   matches $nsb--$nse)\n";
			##### WARNING THIS IS A DISABLE VERSION OF THE PROGRAM
			#$ai[1]=$nqb; $ai[2] =$nqe;
			#$ai[5]=$nsb;  $ai[6] = $nse;
			$t[$i]=(join "\t", @ai);
			splice (@t,$j,1);
			$i--;
			
			print "TABLESIZE, ",scalar(@t),"\n";
			#my $pause=<STDIN>;
			next ILOOP;
			
			
		}
	}
	#my $row=0;
	for (my $i=0; $i < @t ; $i++) {
		$t[$i] =~ s/\r\n/\n/;
		chomp $t[$i];
		my @z=split /\t/, $t[$i];
		my $frow=substr ("000000000$i", -7);
		$z[16]="row$frow" if $i>0;
		#$row++;
		print OUT join ("\t",@z), "\n";	
	}
	close OUT;
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