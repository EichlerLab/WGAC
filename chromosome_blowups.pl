#!/usr/bin/perl
use strict 'vars';
use Getopt::Long;
use File::Copy;

use constant True => 1;
use constant False => 0;
use vars qw(%opt @acc $max_celera $public_length $cbac $sim_auto $sim_sex $sim_min);
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


#build35
my %len;

# TODO: replace hard-coded chromosome lengths with a configuration file.
#-- the GRCh37
%len = qw (chr1 249250621 chr10 135534747 chr11 135006516 chr11_gl000202_random 40103 chr12 133851895 chr13 115169878 chr14 107349540 chr15 102531392 chr16 90354753 chr17 81195210 chr17_ctg5_hap1 1680828 chr17_gl000203_random 37498 chr17_gl000204_random 81310 chr17_gl000205_random 174588 chr17_gl000206_random 41001 chr18 78077248 chr18_gl000207_random 4262 chr19 59128983 chr19_gl000208_random 92689 chr19_gl000209_random 159169 chr1_gl000191_random 106433 chr1_gl000192_random 547496 chr2 243199373 chr20 63025520 chr21 48129895 chr21_gl000210_random 27682 chr22 51304566 chr3 198022430 chr4 191154276 chr4_ctg9_hap1 590426 chr4_gl000193_random 189789 chr4_gl000194_random 191469 chr5 180915260 chr6 171115067 chr6_apd_hap1 4622290 chr6_cox_hap2 4795371 chr6_dbb_hap3 4610396 chr6_mann_hap4 4683263 chr6_mcf_hap5 4833398 chr6_qbl_hap6 4611984 chr6_ssto_hap7 4928567 chr7 159138663 chr7_gl000195_random 182896 chr8 146364022 chr8_gl000196_random 38914 chr8_gl000197_random 37175 chr9 141213431 chr9_gl000198_random 90085 chr9_gl000199_random 169874 chr9_gl000200_random 187035 chr9_gl000201_random 36148 chrM 16571 chrUn_gl000211 166566 chrUn_gl000212 186858 chrUn_gl000213 164239 chrUn_gl000214 137718 chrUn_gl000215 172545 chrUn_gl000216 172294 chrUn_gl000217 172149 chrUn_gl000218 161147 chrUn_gl000219 179198 chrUn_gl000220 161802 chrUn_gl000221 155397 chrUn_gl000222 186861 chrUn_gl000223 180455 chrUn_gl000224 179693 chrUn_gl000225 211173 chrUn_gl000226 15008 chrUn_gl000227 128374 chrUn_gl000228 129120 chrUn_gl000229 19913 chrUn_gl000230 43691 chrUn_gl000231 27386 chrUn_gl000232 40652 chrUn_gl000233 45941 chrUn_gl000234 40531 chrUn_gl000235 34474 chrUn_gl000236 41934 chrUn_gl000237 45867 chrUn_gl000238 39939 chrUn_gl000239 33824 chrUn_gl000240 41933 chrUn_gl000241 42152 chrUn_gl000242 43523 chrUn_gl000243 43341 chrUn_gl000244 39929 chrUn_gl000245 36651 chrUn_gl000246 38154 chrUn_gl000247 36422 chrUn_gl000248 39786 chrUn_gl000249 38502 chrX 155270560 chrY 59373566);

my @chr=qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY);
my @order=qw(chr16 chr19 chr5 chr1 chr2 chr3 chr4 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr17 chr18 chr20 chr21 chr22 chrX chrY);

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

		if(defined($c[6]) && $c[6] ne ''){
		  my $tmp = join "\t", @c;
		  $tmp .= "\n";
		  print OUT $tmp;
		}
		else{
		  print OUT "$c[0]\t$c[1]\t$c[2]\t$c[3]\t0\t9\tgap\n";
		}
	}

	close CHR;
	close OUT;
  
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
		$chrpair{$c[0]}=1;
		$chrpair{$c[4]}=1;
		foreach (0,4) {
			if ($c[$_] eq $chr) {
				$c[$_+1]= int ($c[$_+1]*$scale);
				$c[$_+2]= int ($c[$_+2]*$scale);
				$c[$_+3]=$scaled_length;
			}
		}
		print OUT join ("\t",@c);
	}
	close OUT;
	print "GENERATING LAYOUT FOR ($chr)\n";
	open (OUT , ">chr.layout" ) || die "Can't write chr.layout!\n";
	print OUT "sequence\theader\n";
	foreach my $c (@chr ){
        # Leave blank space for current chr.
		print OUT "$c\t$spot\n" if $c ne $chr && defined $chrpair{$c};
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
}

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
