#!/usr/bin/perl
use strict 'vars';
use Getopt::Long;
use File::Copy;

use constant True => 1;
use constant False => 0;
use vars qw($program $pversion $pdescription $pgenerate);

$program = "$0";
$program =~ s/^.*\///;
### program stats ###
$pversion='1.010727';
$pdescription = "$program (ver:$pversion)  generates parasight views with a single chromosome at center version";
$pgenerate= 'jeff:dnhc genetics:dnh';
### program stats end ###

my %len;

# TODO: replace hard-coded chromosome lengths with a configuration file.
#-- the GRCh37
%len = qw (chr1 249250621 chr10 135534747 chr11 135006516 chr11_gl000202_random 40103 chr12 133851895 chr13 115169878 chr14 107349540 chr15 102531392 chr16 90354753 chr17 81195210 chr17_ctg5_hap1 1680828 chr17_gl000203_random 37498 chr17_gl000204_random 81310 chr17_gl000205_random 174588 chr17_gl000206_random 41001 chr18 78077248 chr18_gl000207_random 4262 chr19 59128983 chr19_gl000208_random 92689 chr19_gl000209_random 159169 chr1_gl000191_random 106433 chr1_gl000192_random 547496 chr2 243199373 chr20 63025520 chr21 48129895 chr21_gl000210_random 27682 chr22 51304566 chr3 198022430 chr4 191154276 chr4_ctg9_hap1 590426 chr4_gl000193_random 189789 chr4_gl000194_random 191469 chr5 180915260 chr6 171115067 chr6_apd_hap1 4622290 chr6_cox_hap2 4795371 chr6_dbb_hap3 4610396 chr6_mann_hap4 4683263 chr6_mcf_hap5 4833398 chr6_qbl_hap6 4611984 chr6_ssto_hap7 4928567 chr7 159138663 chr7_gl000195_random 182896 chr8 146364022 chr8_gl000196_random 38914 chr8_gl000197_random 37175 chr9 141213431 chr9_gl000198_random 90085 chr9_gl000199_random 169874 chr9_gl000200_random 187035 chr9_gl000201_random 36148 chrM 16571 chrUn_gl000211 166566 chrUn_gl000212 186858 chrUn_gl000213 164239 chrUn_gl000214 137718 chrUn_gl000215 172545 chrUn_gl000216 172294 chrUn_gl000217 172149 chrUn_gl000218 161147 chrUn_gl000219 179198 chrUn_gl000220 161802 chrUn_gl000221 155397 chrUn_gl000222 186861 chrUn_gl000223 180455 chrUn_gl000224 179693 chrUn_gl000225 211173 chrUn_gl000226 15008 chrUn_gl000227 128374 chrUn_gl000228 129120 chrUn_gl000229 19913 chrUn_gl000230 43691 chrUn_gl000231 27386 chrUn_gl000232 40652 chrUn_gl000233 45941 chrUn_gl000234 40531 chrUn_gl000235 34474 chrUn_gl000236 41934 chrUn_gl000237 45867 chrUn_gl000238 39939 chrUn_gl000239 33824 chrUn_gl000240 41933 chrUn_gl000241 42152 chrUn_gl000242 43523 chrUn_gl000243 43341 chrUn_gl000244 39929 chrUn_gl000245 36651 chrUn_gl000246 38154 chrUn_gl000247 36422 chrUn_gl000248 39786 chrUn_gl000249 38502 chrX 155270560 chrY 59373566);

my @chr=qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY);
my @order=qw(chr16 chr19 chr5 chr1 chr2 chr3 chr4 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr17 chr18 chr20 chr21 chr22 chrX chrY);

if(scalar(@ARGV) != 6) { die "chromosome_blowups113.pl
    [alignment file]
    [gap/extra file]
    [min bp]
    [min percent identity]
    [chrom to show. all for all chroms, chr1 for chromosome 1 only ......]
    [\"1\" to automate the process (kill the parasight view after generating the ps file. Non \"1\" to keep the window\n";
}
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
