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


## for cow btau4.0
my %len = (
'chr2L' => 23011544,
'chr2LHet' => 368872,
'chr2R' => 21146708,
'chr2RHet' => 3288761,
'chr3L' => 24543557,
'chr3LHet' => 2555491,
'chr3R' => 27905053,
'chr3RHet' => 2517507,
'chr4' => 1351857,
# 'chrM' => 19517,
'chrU' => 10049037,
# 'chrUextra' => 29004656,
'chrX' => 22422827,
'chrXHet' => 204112,
'chrYHet' => 347038           );

my @chr=qw(chr2L chr2LHet chr2R chr2RHet chr3L chr3LHet chr3R chr3RHet chr4 chrU chrX chrXHet chrYHet);
my @order=qw(chr2L chr2LHet chr2R chr2RHet chr3L chr3LHet chr3R chr3RHet chr4 chrU chrX chrXHet chrYHet);

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
my $scaled_length=181800000;
my $screen_width=$scaled_length+3000000;

my $outdir = "starburst.S$min_bp_size.P$min_percent";
mkdir $outdir;
foreach my $chr (@chr) {
	# collect alignments need to draw 
	# just intrachromosomal and interchromosomal for a specific chromosome
        next if(!($chr_show eq "all" || $chr_show eq $chr));
	
	my $scale=0;  # stores the scaling for the chromosomes
	my $spot=1;

print "$len{$chr} and $chr -- \n";	
	$scale= $scaled_length/$len{$chr};
  	print "REMOVING GAPS THAT ARE NOT CHROMOSOME ($chr)\n";
  	open (CHR, $gapfile) || die "Can't read gap file!\n";
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
		print "$c -- \n";
		print OUT "$c\t$spot\n" if $c ne $chr && defined $chrpair{$c};  #leave blank space for current chr
		$spot += $len{$c} + $bp_spacing;
		if ($c eq 'chr3RHet') {
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
print "$chrpair{$c} -- \n";
		if ($c eq $chr) {
			$these_chrs.="$c,$scaled_length:";
			
		} else {
			$these_chrs.="$c,$len{$c}:";
		}
	
	}
	
	my $command= "parasight71.pl  -show $these_chrs -arrange file:chr.layout -extra chr.extras -align chr.alignments";
	
	$command .= " -options '-pair_level=>intra_over_inter, -seq_names_size=>16, -sub_on=>0, -seq_space_paragraph=>120, -seq_space_wrap=>120,-pair_intra_line_on=>1, -pair_inter_line_on=>1, -pair_type_col=>0, -pair_type_col2=>4,-seq_tick_whole=>0, -seq_names_pattern=>chr(\\w+), -seq_tick_on=>0, -seq_tick_b_on=>0, -seq_tick_e_on=>0, -filename_on=>0, -screen_bpwidth =>$screen_width '";
	
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
