#!/usr/bin/perl
use strict 'vars';
use Getopt::Long;
use File::Copy;

use constant True => 1;
use constant False => 0;
use vars qw(%opt @acc $max_celera $public_length $cbac $sim_auto $sim_sex $sim_min);

###MODIFICATIONS###
### 01-07-27  altered the program to remove temp files when run in non-batch mode###
### 2003-04-02 altered from Jeff's chromosome_blowups11.pl to fit the needs of Xinwei's chr19 project 
###2003-6-18 modified to add more extras on the track: they should have at least 7 fields and added in the file Npositions_chr19.tbl, in which gap extras only have 6 fields.
###2003/12/19 modified to draw chr2/4
##	2007-6-20 to modify it for the pantro2, to add the %len


use vars qw($program $pversion $pdescription $pgenerate);
$program = "$0";
$program =~ s/^.*\///;
### program stats ###
$pversion='1.010727';
$pdescription = "$program (ver:$pversion)  generates parasight views with a single chromosome at center version";
$pgenerate= 'jeff:dnhc genetics:dnh';
### program stats end ###
#print "usage: $program -in [path] -out [path] [options]\n";
#

#oo33
my %len = qw (chr1 287033928 chr1_random 930219 chr2 253614609 chr3 227588907 chr3_random 418273 chr4 209899732 chr4_random 635801 chr5 201269662 chr5_random 1140379 chr6 192062348 chr7 172546273 chr7_random 172546273 chr8 164220812 chr8_random 2266 chr9 152859692 chr9_random 987849 chrX 159673156 chrX_random 3031746 chrY 50860226 chr10 146852389 chr10_random 2071655 chr11 154564754 chr11_random 400000 chr12 154076372 chr12_random 154076372 chr13 118740136 chr13_random 118740136 chr14 106660743 chr14_random 200000 chr15 108032665 chr15_random 274431 chr16 108615791 chr17 89876491 chr17_random 707248 chr18 90206333 chr18_random 358746 chr19 72301379 chr19_random 478517 chr20 64863847 chr21 44907570 chr22 47748585 chrNA_random 54262405 chrUL_random 30253343); 

#build29 with publish quality chr7
%len = qw (chr1 245933730 chr10 138509279 chr10_random 1347815 chr11 134078264 chr11_random 4496094 chr12 134180709 chr12_random 1886052 chr13 113614620 chr13_random 1867447 chr14 105387637 chr14_random 752709 chr15 100020159 chr15_random 3566954 chr16 93247745 chr16_random 1612449 chr17 81112292 chr17_random 4426481 chr18 80862210 chr18_random 4964069 chr19 73000771 chr19_random 464170 chr1_random 2554216 chr2 246757188 chr2_random 4622209 chr3 204706827 chr3_random 7996578 chr4 196955444 chr4_random 5809204 chr5 189570559 chr5_random 4123556 chr6 182825515 chr6_random 5311953 chr7 157432593 chr8 148433236 chr8_random 4105857 chr9 129119237 chr9_random 785660 chrUn_random 12510748 chrX 151820622 chrX_random 2012070 chrY 58368225 chr20 62842997 chr21 44626493 chr22 47748585) ;

#build31 with publish quality chr19
%len = qw (chr1 244258774 chr10 134647902 chr10_random 141075 chr11 136521022 chr11_random 1362657 chr12 133382389 chr12_random 1258214 chr13 111298136 chr14 101218245 chr14_random 64039 chr15 96598362 chr15_random 817318 chr16 91211881 chr16_random 382765 chr17 84296999 chr17_random 1434149 chr18 78067305 chr19 58789685 chr1_random 584119 chr2 241996787 chr20 62802940 chr21 44709625 chr22 47748585 chr2_random 160669 chr3 199558344 chr3_random 637829 chr4 191669278 chr4_random 953437 chr5 181762559 chr6 170670676 chr6_random 9046941 chr7 157432593 chr8 146305119 chr8_random 526343 chr9 132877114 chr9_random 36248 chrUn_random 1774146 chrX 151567156 chrY 50360226);

#build34 with publish quality chr2/4
%len = qw (chr1 246127941 chr2 243615958 chr3 199344050 chr4 191731959 chr5 181034922 chr6 170914576 chr7 158545518 chr8 146308819 chr9 136372045 chrX 153692391 chrY 50286555 chr1_random 6515988 chr2_random 1104831 chr3_random 749256 chr4_random 648024 chr5_random 143687 chr6_random 2055751 chr7_random 632637 chr8_random 1499381 chr9_random 2766341 chrX_random 3403558 chr10 135037215 chr11 134482954 chr12 132078379 chr13 113042980 chr14 105311216 chr15 100256656 chr16 90041932 chr17 81860266 chr18 76115139 chr19 63811651 chr20 63741868 chr21 46976097 chr22 49396972 chr10_random 1043775 chr13_random 189598 chr15_random 1132826 chr17_random 2549222 chr18_random 4262 chr19_random 92689 chrUn_random 3349625);

#build35
%len = qw (chr1 245522847 chr1_random 3897131 chr2 243018229 chr2_random 418158 chr3 199505740 chr3_random 970716 chr4 191411218 chr4_random 1030282 chr5 180857866 chr5_random 143687 chr6 170975699 chr6_random 1875562 chr7 158628139 chr7_random 778964 chr8 146274826 chr8_random 943810 chr9 138429268 chr9_random 1312665 chr10 135413628 chr10_random 113275 chr11 134452384 chr12 132449811 chr12_random 466818 chr13 114142980 chr13_random 186858 chr14 106368585 chr15 100338915 chr15_random 784346 chr16 88827254 chr16_random 105485 chr17 78774742 chr17_random 2618010 chr18 76117153 chr18_random 4262 chr19 63811651 chr19_random 301858 chr20 62435964 chr21 46944323 chr22 49554710 chr22_random 257318 chrX 154824264 chrX_random 1719168 chrY 57701691);

## pantro2
 %len = (
'chr1' => 229974691,
'chr10' => 135001995,
'chr10_random' => 8402541,
'chr11' => 134204764,
'chr11_random' => 8412303,
'chr12' => 135371336,
'chr12_random' => 2259969,
'chr13' => 115868456,
'chr13_random' => 9231045,
'chr14' => 107349158,
'chr14_random' => 2108736,
'chr15' => 100063422,
'chr15_random' => 3087076,
'chr16' => 90682376,
'chr16_random' => 6370548,
'chr17' => 83384210,
'chr17_random' => 5078517,
'chr18' => 77261746,
'chr18_random' => 1841920,
'chr19' => 64473437,
'chr19_random' => 2407237,
'chr1_random' => 9420409,
'chr20' => 62293572,
'chr20_random' => 1792361,
'chr21' => 46489110,
'chr22' => 50165558,
'chr22_random' => 1182457,
'chr2a' => 114460064,
'chr2a_random' => 3052259,
'chr2b' => 248603653,
'chr2b_random' => 2186977,
'chr3' => 203962478,
'chr3_random' => 3517036,
'chr4' => 194897272,
'chr4_random' => 6711082,
'chr5' => 183994906,
'chr5_random' => 3159943,
'chr6' => 173908612,
'chr6_random' => 9388360,
'chr7' => 160261443,
'chr7_random' => 7240870,
'chr8' => 145085868,
'chr8_random' => 7291619,
'chr9' => 138509991,
'chr9_random' => 7733331,
'chrUn' => 58616431,
'chrX' => 155361357,
'chrX_random' => 3548706,
'chrY' => 23952694,
'chrY_random' => 772887,
'Y1'		=>  115015,
'Y2' => 235662,
'Y3' => 66393,
'Y4' => 554625
                        );

my @chr=qw(chr1 chr2a chr2b chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY Y1 Y2 Y3 Y4 chr1_random chr2a_random chr2b_random chr3_random chr4_random chr5_random chr6_random chr7_random chr8_random chr9_random chr10_random chr11_random chr12_random chr13_random chr14_random chr15_random chr16_random chr17_random chr18_random chr19_random chr20_random chr21_random chr22_random chrX_random chrY_random);
my @order=qw(chr16 chr19 chr5 chr1 chr2a chr2b chr3 chr4 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr17 chr18 chr20 chr21 chr22 chrX chrY Y1 Y2 Y3 Y4  chr1_random chr2a_random chr2b_random chr3_random chr4_random chr5_random chr6_random chr7_random chr8_random chr9_random chr10_random chr11_random chr12_random chr13_random chr14_random chr15_random chr16_random chr17_random chr18_random chr19_random chr20_random chr21_random chr22_random chrX_random chrY_random);

if(scalar(@ARGV) != 6) { die "chromosome_blowups113.pl \n\t[alignment file]\n\t[Gap/extra fiel]\n\t[min bp]\n\t[min pcent identity]\n\t[chrom to show. all for all chroms, chr1 for chromosome 1 only ......]\n\t[\"1\" to automate the process (kill the parasight view after generating the ps file. Non \"1\" to keep the window\n";}
my $alignfile = $ARGV[0];
my $gapfile = $ARGV[1];
my $chr_show = $ARGV[4];

my $min_bp_size=$ARGV[2];
my $min_duppos_size=10000;
die "usage: chromsome_blowups1.pl [min bases] [min frac identity]\n" if !$ARGV[0] || !$ARGV[1];
my $min_percent=$ARGV[3];
my $automate = $ARGV[5];

my $bp_spacing=10000000;
my $scaled_length=1818000000;
my $screen_width=$scaled_length+3000000;

my $outdir = "starburst.S$min_bp_size.P$min_percent";
mkdir $outdir;
foreach my $chr (@chr) {
	# collect alignments need to draw 
	# just intrachromosomal and interchromosomal for a specific chromosome
        next if(!($chr_show eq "all" || $chr_show eq $chr));
	
	my $scale=0;  # stores the scaling for the chromosomes
	my $spot=1;
	
	$scale= $scaled_length/$len{$chr};
  	print "REMOVING GAPS THAT ARE NOT CHROMOSOME ($chr)\n";
  	open (CHR, $gapfile) || die "Can't read Npositions_chr19.tbl!\n";
  	open (OUT, ">chr.extras") || die "Can't write chr.extras!\n";
 	my $head=<CHR>;
 	print OUT "seq\tbegin\tend\tcolor\toffset\twidth\ttype\n";
  	my %chrpair=();
  	while (<CHR>) {
  		s/\r\n/\n/;
		chomp;
		my @c = split /\t/;
		next if $c[0] ne $chr && $c[3] eq 'white';  
		if ($c[0] eq $chr) {
			$c[1]= int ($c[1]*$scale);
			$c[2]= int ($c[2]*$scale);
		
		}

		#####the below was added on 6/18/2003
		if(defined($c[6]) && $c[6] ne ''){
		  my $tmp = join "\t", @c;
		  $tmp .= "\n";
		  print OUT $tmp;
		}
		else{
		##### the above was added on 6/18/2003

		  print OUT "$c[0]\t$c[1]\t$c[2]\t$c[3]\t0\t9\tgap\n";
		}
		#white means gaps (centromeres and parms are purple)
		#for (my $i=0; $i< @c; $i++) {
		#	print "$i)$c[$i] ";
		#}
		#print "\n";
		#my $pause=<STDIN>;
	}
	close CHR;
	

##### Comment out hotspots. May use later ######
#  	open (CHR, "oo.2N.weld10kb.join.all.cull.hotspots.reformat.span.nooverlap") || die "Can't read oo.2N.weld10kb.join.all.cull.hotspots.reformat.span.nooverla!\n";
#	my $head=<CHR>;
# 	my %chrpair=();
# 	while (<CHR>) {
# 		s/\r\n/\n/;
#		chomp;
#		my @c = split /\t/;
#		next if $c[0] ne $chr && $c[3] eq 'white';  
#		if ($c[0] eq $chr) {
#			$c[1]= int ($c[1]*$scale);
#			$c[2]= int ($c[2]*$scale);
#			$c[3]='orange';
#			$c[4]='red';
#		}
#		print OUT "$c[0]\t$c[1]\t$c[2]\torange\t10\t9\thotspot\n";
#	}
	close OUT;
#	close CHR;
	
#  die "";
  #print "COLLECTING CELERA DUP POSITIVE REGIONS TO DISPLAY ($chr)\n";
 # open (CHR, "DupPosNonredundant.txt.nooverlap") || die "Can't read DupPosNonredundant.txt.nooverlap!\n";
 # open (OUT, ">chr.duppos.extras") || die "Can't writechr.duppos.extras\n";
 # my $head =<CHR>;
 # print OUT "seq\tbegin\tend\tcolor\toffset\twidth\n";
 # while (<CHR>) {
#		my @c = split /\t/;
#
#		next if $c[0] ne $chr ;
#		next if $c[2]-$c[1] < $min_duppos_size;
#		if ($c[0] eq $chr) {
#			$c[1]= int ($c[1]*$scale);
#			$c[2]= int ($c[2]*$scale);
#		
#		}
#		#for (my $i=0; $i< @c; $i++) {
		#	print "$i)$c[$i] ";
		#}
		#print "\n";
#		
#		print OUT "$c[0]\t$c[1]\t$c[2]\tblack\t15\t12\n";
#		#my $pause=<STDIN>;
#	}
#	close OUT;
  
  print "REMOVING ALIGNMENTS THAT ARE NOT CHROMOSOME ($chr)\n";
  open (CHR, $alignfile) || die "Can't read alignment.pieces!\n";
  open (OUT, ">chr.alignments") || die "Can't write chr.alignments!\n";
  my $head =<CHR>;
  print OUT $head;
  while (<CHR>) {
		my @c = split /\t/;
		next if $c[27] < $min_percent;
		next if $c[0] ne $chr && $c[4] ne $chr;
		next if $c[22] < $min_bp_size;
		next if $c[0] =~ /random/;
		next if $c[4] =~ /random/;
	 next if $c[0] =~ /chrUn/;
		 next if $c[4] =~ /chrUn/;
		$chrpair{$c[0]}=1;
		$chrpair{$c[4]}=1;
		foreach (0,4) {
			if ($c[$_] eq $chr) {
				$c[$_+1]= int ($c[$_+1]*$scale);
				$c[$_+2]= int ($c[$_+2]*$scale);
				$c[$_+3]=$scaled_length;
			}
		}
		#for (my $i=0; $i< @c; $i++) {
		#	print "$i)$c[$i] ";
		#}
		#print "\n";
		print OUT join ("\t",@c);
		#my $pause=<STDIN>;
	}
	close OUT;
	print "GENERATING LAYOUT FOR ($chr)\n";
	open (OUT , ">chr.layout" ) || die "Can't write chr.layout!\n";
	print OUT "sequence\theader\n";
	foreach my $c (@chr ){
		print OUT "$c\t$spot\n" if $c ne $chr && defined $chrpair{$c};  #leave blank space for current chr
		$spot += $len{$c} + $bp_spacing;
		if ($c eq 'chr8') {
			#inset the middle line#
			print OUT "NEWLINE\n";
			print OUT "$chr\t1\n";
			print OUT "NEWLINE\n";
			$spot=1;
			
		}

	}
	print "GENERATING SHOW LINE FOR ($chr)\n";
	my $these_chrs='';
	foreach my $c (@chr) {
		next if !defined $chrpair{$c};
		if ($c eq $chr) {
			$these_chrs.="$c,$scaled_length:";
			
		} else {
			$these_chrs.="$c,$len{$c}:";
		}
	
	}
	
	my $command= "parasight71.pl  -show $these_chrs -arrange file:chr.layout -extra chr.extras -align chr.alignments";
	
	$command .= " -options '-pair_level=>intra_over_inter, -seq_names_size=>16, -sub_on=>0, -seq_space_paragraph=>120, -seq_space_wrap=>120,-pair_intra_line_on=>1, -pair_inter_line_on=>1, -pair_type_col=>0, -pair_type_col2=>4,-seq_tick_whole=>0, -seq_names_pattern=>chr(\\w+), -seq_tick_on=>0,-seq_tick_b_on=>0, -seq_tick_e_on=>0, -filename_on=>0, -screen_bpwidth =>$screen_width '";
	
	$command .= " -precode '\$canvas->raise(inter); \$canvas->raise(intra); \$canvas->raise(seqn); &print_screen(0); '  ";
	
	$command .=" -die" if($automate == 1);
	
	print "$command\n";
	system $command;
	my $kb=int($min_bp_size/1000);
	my $perc=int($min_percent*1000)/10;
	system "mv screen.ps $outdir/$chr"."_$kb"."kb_$perc"."perc.ps";
	
	
#		." -options '-seq_tick_whole=>0,-extra_label_size=>12,-filename_on=>0,-text_text=>$name, 
#-text_offset=>250, -graph1_label_size=12,-text_offset_h=>280,-text_color=>green,-graph2_label_on=>0,-graph1_label_decimal=>3,-gscale_indent=>-60,-graph1_max=>1,-gscale_on=>1,-graph1_min=>0.98,-screen_indent_l=>80,-screen_bpwidth=>$width,-extra_label_on=>1,-seq_tick_bp=>10000,-extra_label_col=>$label_column, -seq_names_offset=>7, -seq_names_size=>13,-seq_names_color=>darkgrey, -extra_label_color=>black,-window_width=>700, -window_height=>400, -pair_intra_line_on=>1,-sub_on=>0,-seq_tick_e_label_offset=>-4,-seq_tick_e_length=>1,-seq_tick_e_label_offset_h=>12, -seq_tick_e_label_size=>11'";




#	my $label_column=8;
#	my $width=$tot_seq;
#	my $command= "parasight71.pl  -show $opt{'out'}.show -graph1 $opt{'out'}.graph1 -arrange file:$opt{'out'}.arrange -extra $opt{'out'}.extra"
#		." -options '-seq_tick_whole=>0,-extra_label_size=>12,-filename_on=>0,-text_text=>$name, -text_offset=>250, -graph1_label_size=12,-text_offset_h=>280,-text_color=>green,-graph2_label_on=>0,-graph1_label_decimal=>3,-gscale_indent=>-60,-graph1_max=>1,-gscale_on=>1,-graph1_min=>0.98,-screen_indent_l=>80,-screen_bpwidth=>$width,-extra_label_on=>1,-seq_tick_bp=>10000,-extra_label_col=>$label_column, -seq_names_offset=>7, -seq_names_size=>13,-seq_names_color=>darkgrey, -extra_label_color=>black,-window_width=>700, -window_height=>400, -pair_intra_line_on=>1,-sub_on=>0,-seq_tick_e_label_offset=>-4,-seq_tick_e_length=>1,-seq_tick_e_label_offset_h=>12, -seq_tick_e_label_size=>11'";
#	if ($opt{'psonly'} ) {
#		$command .= " -precode '&print_screen(0);  ' -die ";
#	} else {
#		#$command .= " -precode 'unlink \"$opt{out}.align\"; unlink \"$opt{out}.extras\";' " ;
#	}
#	print "$command\n";
#	system $command;
#	#die "";


}
###################################
####subroutines####################
	
sub lookup {
	my $input=shift;
	my $acc=$input;
	$acc=~s/\.\d+$//;
	my $last2=substr($acc,-2);
	print "$last2\n";
	my $dir="$opt{'db'}/$last2";
	opendir(DIR,$dir) || die "Can't read directory ($dir)!\n";
	my @accdir=grep { /$input/} readdir DIR;
	close DIR;
	print "No directory found matching ($input)!\n" if @accdir==0;
	print  "Too many matching ($input) => ",join(" ",@accdir), "\n"  if @accdir >1;
	$dir .= "/$accdir[0]";
	return ($dir,$accdir[0]) if @accdir==1; 
	#my $pause=<STDIN>;
	next ACC_LOOP;  #accdir[0] really version#
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


=head2


































































































































	$max_celera= &max_position_celera($dir);
	print "MAX_POS_CELERA:$max_celera\n";
	&generate_pairwise( $dir);
	

	open (FASTA,"$dir/$ver") || die "Can't open fasta ($dir/$ver)!\n";
	my @ns=&fasta_ns_positions (-filehandle=> \*FASTA,-min_size=>20);
	foreach (@ns) {
		#$$_{'end'};
		print EXTRA "$ver\t$$_{'begin'}\t$$_{'end'}\tN \tblack\t-7\t5\t\t\n";
	}
	close FASTA;
	
	
	
	###GET MAX_#READS IN WINDOW AND CALCULATE SCALE###
	my $max_win=0;
	open (WIN, "$dir/$ver.win") || die "Can't read file ($dir/$ver.win)!\n";
	while (<WIN>) {
		chomp;
		s/ScaCtg/C/;
		my @c=split / +/;
		my $contig=shift @c;
		for (my $i=0; $i<@c; $i++) { 
			$max_win=$c[$i] if $c[$i] > $max_win;
		}
	}
	close WIN;
	#####set up certain things#####
	my $max=$max_win;
	$max=105 if $max <100;
	my $scale= 200/$max;
	my $ypos= 10 + $opt{'x_thresh'} * $scale;
	my $lname='SEX';
	$lname=$opt{'type'} if $opt{'type'} ne 'BOTH';
	print EXTRA "$cbac\t1\t$max_celera\t$lname:$opt{'x_thresh'} \tpink\t$ypos\t1\tX:$opt{'x_thresh'} \n" if $opt{'type'} =~ /BOTH|SEX|X|Y/;
	$ypos= 10 + $opt{'a_thresh'} * $scale;
	$lname='AUTO';
	$lname=$opt{'type'} if $opt{'type'} ne 'BOTH';
	print EXTRA "$cbac\t1\t$max_celera\t$lname:$opt{'a_thresh'} \tred\t$ypos\t1\tAUTO:$opt{'a_thresh'} \n" if $opt{'type'} !~ /SEX|X|Y/;
	my $y=100;
	while ($y < $max) {
		$ypos= 10 + $y* $scale;
		print EXTRA "$cbac\t1\t$max_celera\t$y \tblack\t$ypos\t1\t$y \n";
		$y*=2;
	}
	$ypos= 10 + $max_win * $scale;
	print EXTRA "$cbac\t1\t$max_celera\tMAX:$max \tdarkgrey\t$ypos\t1\tMAX:$max_win \n";
	
	open (WIN, "$dir/$ver.win") || die "Can't read file ($dir/$ver.win)!\n";
	while (<WIN>) {
		chomp;
		my @c=split / +/;
		my $contig=shift @c;
		$contig =~s/ScaCtg/C/;
		$contig =~s/://;
		$contig =~s/^>//;
		#print "$contig\n";
		#print join ("X", keys %ctg), "\n";
		for (my $i=0; $i<@c; $i++) {
			my $begin= $i*1000+$ctg{$contig}+1;
			my $end=$begin+4999;
			my $yposition= 10 + $c[$i]* $scale;
			my $color='grey';
			if ($c[$i]>$opt{'x_thresh'}) {
				if ($opt{'type'} eq 'BOTH') {
					$color='pink' ;
				} elsif ($opt{'type'} =~/SEX|X|Y/) { 
					$color='red';
				}
			}
			$color='red' if $c[$i]>$opt{'a_thresh'};
			print EXTRA "$cbac\t$begin\t$end\t$c[$i]\t$color\t$yposition\t2\t\t\n";
		}
	}
	
	my $width=100000;
	$width=$max_celera+1000 if $max_celera+1000>$width;
	$width=$public_length+1000 if $public_length+1000>$width;
	my %colorhash = ('SINE/Alu' => 'blue', 'SINE/MIR' => 'blue', 'LINE/L1' => 'medium orchid',
			'LINE/L2' => 'medium orchid', 'DNA/MER1_type' => 'light salmon',
			'DNA/MER2_type' => 'light salmon', 'LTR/ERVL' => 'cyan', 'LTR/MaLR' => 'cyan',
			'Low_complexity' => 'red', 'Simple_repeat' => 'red', 'Satellite/telo' => 'medium orchid',
			'LTR/ERV1' => 'cyan', 'LINE/CR1' => 'medium orchid', 'snRNA' => 'orange', 'LTR/ERV' => 'cyan',
			'LTR/MER21-group' => 'cyan'
	);

	my %offsethash = ('SINE/Alu' => '-12', 'SINE/MIR' => '-12', 'LINE/L1' => '-17',
			'LINE/L2' => '-17', 'DNA/MER1_type' => '-26',
			'DNA/MER2_type' => '-26', 'LTR/ERVL' => '-26', 'LTR/MaLR' => '-26',
			'Low_complexity' => '-21', 'Simple_repeat' => '-21', 'Satellite/telo' => '-26',
			'LTR/ERV1' => '-26', 'LINE/CR1' => '-17', 'snRNA' => '-26', 'LTR/ERV' => '-26',
			'LTR/MER21-group' => '-26'
	);
	
	print EXTRA "$ver\t1\t1\tSINES \twhite\t-12\t4\tSINES\t\t\t\n";
	print EXTRA "$ver\t1\t1\tLINES \twhite\t-17\t4\tLINES\t\t\t\n";
	print EXTRA "$ver\t1\t1\tSIMPLE/LOW \twhite\t-21\t4\tSIMPLE/LOW\t\t\t\n";
	print EXTRA "$ver\t1\t1\tOTHER \twhite\t-26\t4\tOTHER\t\t$\t\n";
	print EXTRA "$ver\t1\t1\t<3% \twhite\t4\t4\t<3%t\t\t\n";

	open (RMOUT, "$dir/$ver.out") || die "Can't read file ($dir/$ver.out)!\n";
	while (<RMOUT>) {
      s/\r\n/\n/;
      chomp;
      s/^ +//;
      next if !(/^\d/);
      my @c = split / +/;
		my $color='green';
		my $offset='-9';
		my $orientation='';
		$color=$colorhash{$c[10]} if defined ($colorhash{$c[10]});
		$offset=$offsethash{$c[10]} if defined ($offsethash{$c[10]});

	#RepeatMasker column 9 (orientation) conversion to "F/R" from "+/C"
		if ($c[8] eq "+") {
			$orientation = "F";
			$c[19] = $orientation;
			}
		elsif ($c[8] eq "C"){
			$orientation = "R";
			$c[19] = $orientation;
			}
		print EXTRA "$ver\t$c[5]\t$c[6]\t\t$color\t$offset\t4\t\t$c[9]:$c[10]\t$c[1]\t\n";
		if ( $c[1] < 3) {
			next if $c[10] =~/Simple|complexity/;
			print EXTRA "$ver\t$c[5]\t$c[6]\t\thotpink\t7\t4\t\t$c[9]:$c[10]\t$c[1]\t\n";
		}
		
	}
	close RMOUT;
	close EXTRA;
	
	#die "";
	
	my $label_column=3;
	$label_column=7 if $opt{'label'};
	my $command= "parasight71.pl -align $opt{'out'}.align -extra $opt{'out'}.extras -show $ver:$cbac "
	." -options '-screen_indent_l=>80,-screen_bpwidth=>$width,-extra_label_on=>1,-seq_tick_bp=>10000,-extra_label_col=>$label_column, -seq_names_offset=>-14, -seq_names_size=>26,-seq_names_color=>darkgrey, -extra_label_color=>black,-window_width=>700, -window_height=>400, -pair_level=>intra_over_inter, -pair_intra_line_on=>1,-sub_on=>0,-seq_tick_e_label_offset=>-4,-seq_tick_e_length=>1,-seq_tick_e_label_offset_h=>12, -seq_tick_e_label_size=>11'";
if ($opt{'psonly'} ) {
	$command .= " -precode ' &print_screen(0)' -die ";
} else {
	$command .= " -precode 'unlink \"$opt{out}.align\"; unlink \"$opt{out}.extras\";' " ;
}
	print "$command\n";
	system $command;

}

#$opt{'seq_tick_e_on'}=$true if !defined $opt{'seq_tick_e_on'};
#$opt{'seq_tick_e_offset'}=0 if !defined $opt{'seq_tick_e_offset'};
#$opt{'seq_tick_e_width'}=2 if !defined $opt{'seq_tick_e_width'};
#$opt{'seq_tick_e_length'}=10 if !defined $opt{'seq_tick_e_length'};
#$opt{'seq_tick_e_color'}='black' if !defined $opt{'seq_tick_e_color'};
#$opt{'seq_tick_e_label_on'}=$true if !defined $opt{'seq_tick_e_label_on'};
#$opt{'seq_tick_e_label_multiplier'}=0.001 if !defined $opt{'seq_tick_e_label_multiplier'};
#$opt{'seq_tick_e_label_size'}=9 if !defined $opt{'seq_tick_e_label_size'};
#$opt{'seq_tick_e_label_color'}='black' if !defined $opt{'seq_tick_e_label_color'};
#$opt{'seq_tick_e_label_offset'}=2 if !defined $opt{'seq_tick_e_label_offset'};
#3$opt{'seq_tick_e_label_offset_h'}=0 if !defined $opt{'seq_tick_e_label_offset_h'};


sub generate_pairwise {
	my $dir = shift;
	my $ver = $dir;
	$ver =~ s/^.*\///;
	print "VER:$ver\n";
	
	my @colors=("darkgreen","darkblue","purple","blue","darkred");
	open (ALIGN, ">$opt{'out'}.align") || die "Can't open alignment file ($opt{'out'}.align)!\n";
	print ALIGN "name1\tbegin1\tend1\tlength1\tname2\tbegin2\tend2\tlength2\tcolor\n";
	open (B2P, "$dir/$ver.b2p") || die "Can't read file ($dir/$ver.b2p)!\n";
	while (<B2P>) {
		chomp;
		my @c=split / +/;
		$c[0]= "ubac$c[0]";
		$cbac=$c[0];
		$cbac=~s/_.*$//;
		$c[1]=$ver;
		$c[2]+=1;
		$c[4]+=1;
		my $orient=pop @c;
		print "ORIENT:($orient)\n";
		if ($orient eq 'R') {
			print "FLIPPING\n";
			($c[4],$c[5]) = ($c[5],$c[4]);
		}
		print join(":",@c), "\n";
		push @colors, shift @colors;
		print ALIGN "$cbac\t$c[2]\t$c[3]\t$max_celera\t$c[1]\t$c[4]\t$c[5]\t$public_length\t$colors[0]\tdummy\n";
	}
	close ALIGN;
}

sub fasta_ns_positions {
	my %args=(-filehandle=>\*STDIN, -min_size=>1 ,@_);
	my ($true,$false)=(1,0);
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

sub max_position_celera {
	my $dir = shift;
	my $ver = $dir;
	my $max_length=0;
	$ver =~ s/^.*\///;
	print "VER:$ver\n";
	open (B2P, "$dir/$ver.b2p") || die "Can't read file ($dir/$ver.b2p)!\n";
	while (<B2P>) {
		chomp;
		print;
		my @c=split / +/;
		$max_length=$c[3] if $c[3]>$max_length;
	}
	close B2P;
	open (B2C, "$dir/$ver.b2c") || die "Can't read file ($dir/$ver.b2c)!\n";
	while (<B2C>) {
		chomp;
		my @c=split / +/;
		$max_length=$c[2] if $c[2]>$max_length;
	}
	close B2C;
	return $max_length;
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

Copyright (C) 2001 Jeff Bailey. Extremely Large Monetary Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut
