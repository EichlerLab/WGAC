#!/bin/env perl
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

if(scalar(@ARGV) != 7) {
    die "chromosome_blowups.pl
    [alignment file]
    [gap/extra file]
    [min bp]
    [min percent identity]
    [chrom to show. all for all chroms, chr1 for chromosome 1 only ......]
    [\"1\" to automate the process (kill the parasight view after generating the ps file. Non \"1\" to keep the window
    [chromosome length and order file]
";
}
my $alignfile = $ARGV[0];
my $gapfile = $ARGV[1];
my $chr_show = $ARGV[4];

my $min_bp_size=$ARGV[2];
my $min_duppos_size=10000;
my $min_percent=$ARGV[3];
my $automate = $ARGV[5];
my $chromosome_lengths_file = $ARGV[6];

my $bp_spacing=150000000;
my $scaled_length=1818000000;
my $screen_width=$scaled_length+3000000;
# my $bp_spacing=20000000;
# my $scaled_length=1820000000;
# my $screen_width=$scaled_length+30000000;

# Read chromosome lengths keeping track of lengths by chromosome name and the
# order of the chromosomes themselves.
my %len = ();
my @chr= ();
my $chromosome;
my $length;
open(CHR, $chromosome_lengths_file) ||
    die "Can't read chromosome lengths file: $chromosome_lengths_file\n";
while (<CHR>) {
    chomp;
    ($chromosome, $length) = split(/\t/);

    # Save the length by chromosome name.
    $len{$chromosome} = $length;

    # Save the chromosome name in the order it is read in.
    push(@chr, $chromosome);
}
close(CHR);

# Define indices for necessary columns in the input file.
my $min_percent_index = 10;
my $min_bp_size_index = 8;
#my $min_percent_index = 27;
#my $min_bp_size_index = 22;
my $chr_first_index = 0;
my $chr_second_index = 4;

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
        next if $c[$min_percent_index] < $min_percent;
        next if $c[$chr_first_index] ne $chr && $c[$chr_second_index] ne $chr;
        next if $c[$min_bp_size_index] < $min_bp_size;
        next if $c[$chr_first_index] =~ /random/;
        next if $c[$chr_second_index] =~ /random/;
        $chrpair{$c[$chr_first_index]}=1;
        $chrpair{$c[$chr_second_index]}=1;
        foreach ($chr_first_index,$chr_second_index) {
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

    if (length($these_chrs) == 0) {
        print "Couldn't find any data for show line. Skipping $chr.\n";
        next;
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
