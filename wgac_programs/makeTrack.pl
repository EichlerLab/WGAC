#!/usr/bin/env perl
# 
# 1st arg:  output of previous step in WGAC pipeline (oo.weld10kb.join.all.cull)
# output:  a new table/track on the UCSC Genome Browser

use strict;

my ($infile) = @ARGV;

open (INFILE, "$infile") || die "cannot read input file $infile:  $!";
$_ = <INFILE>;
chomp;
my @allColumnNames = split /\t/;

my %colname2Pos; # hashtable translates names of columns to positions (0,1,2...)
my $index = 0;
%colname2Pos = map {$_ => $index++} @allColumnNames;

@allColumnNames = (); # don't need this array anymore

unless (   exists $colname2Pos{'QNAME'}
        && exists $colname2Pos{'QB'}
        && exists $colname2Pos{'QE'}
        && exists $colname2Pos{'QLEN'}
        && exists $colname2Pos{'SNAME'}
        && exists $colname2Pos{'SB'}
        && exists $colname2Pos{'SE'}
        && exists $colname2Pos{'SLEN'}
		&& exists $colname2Pos{'FILE'} # to be stored on web server
        && exists $colname2Pos{'END'} # length of alignment (incl. gaps & Ns) = alignL 
        && exists $colname2Pos{'indel_N'}
        && exists $colname2Pos{'indel_S'}
        && exists $colname2Pos{'base_S'} # num of alignment positions between real bases (not Ns, not gaps), i.e. matches plus mismatches = alignB
        && exists $colname2Pos{'base_Match'}
        && exists $colname2Pos{'base_Mis'}
        && exists $colname2Pos{'transitions'}
        && exists $colname2Pos{'transversions'}
        && exists $colname2Pos{'per_sim'} # fracMatch
        && exists $colname2Pos{'per_sim_indel'} # fracMatchIndel
        && exists $colname2Pos{'K_jc'}
        && exists $colname2Pos{'k_kimura'} ) {
	die "input file has one or more missing columns \n";		
}

my %rowsPrinted; 
my $rowNumber = 0;

while (<INFILE>) {
	next unless /\S/;
	chomp;
	my @columnvalues = split /\t/;
	
	# skip if zero length
	next unless $columnvalues[$colname2Pos{'QLEN'}];
	next unless $columnvalues[$colname2Pos{'SLEN'}];

	
	my $chrom = $columnvalues[$colname2Pos{'QNAME'}];
	my $chromStart = $columnvalues[$colname2Pos{'QB'}];
	my $chromEnd = $columnvalues[$colname2Pos{'QE'}];
	
	my $otherChrom = $columnvalues[$colname2Pos{'SNAME'}];
	my $otherStart;
	my $otherEnd;
	my $strand;
	if ($columnvalues[$colname2Pos{'SB'}] < $columnvalues[$colname2Pos{'SE'}]) {
		$otherStart = $columnvalues[$colname2Pos{'SB'}];
		$otherEnd = $columnvalues[$colname2Pos{'SE'}];
		$strand = '+';
	} else {
		$otherStart = $columnvalues[$colname2Pos{'SE'}];
		$otherEnd = $columnvalues[$colname2Pos{'SB'}];
		$strand = '_';
	}
	
	# convert to UCSC convention
	--$chromStart;
	--$otherStart;
	
	my $rowID = "$chrom:$chromStart:$chromEnd:$otherChrom:$otherStart:$otherEnd";
	unless (exists $rowsPrinted{$rowID}) {
		++$rowsPrinted{$rowID}; # remember that we've seen this
		++$rowNumber;
		printRow($rowID, $strand, $rowNumber, \@columnvalues);
	}
	
	# also print the inverse (using the same rowNumber), unless it would be a duplicate row
	$rowID = "$otherChrom:$otherStart:$otherEnd:$chrom:$chromStart:$chromEnd";
	unless (exists $rowsPrinted{$rowID}) {
		++$rowsPrinted{$rowID};
		printRow($rowID, $strand, $rowNumber, \@columnvalues);
	}
	
	
}
close INFILE;
exit;


# print 
# chrom, start, end, otherChrom, start, end, strand, otherChrom:start, id number, alignfile
# alignL = length of alignment (incl. gaps & Ns)
# 

sub printRow {
	my $rowID = shift @_;
	my $strand = shift @_;
	my $rowNumber = shift @_;
	my $refcolumnValues = shift @_;
	
	my ($chrom, $chromStart, $chromEnd, $otherChrom, $otherStart, $otherEnd) = split (':', $rowID);
	
	print "$chrom\t$chromStart\t$chromEnd\t";
	print "$otherChrom:$otherStart\t"; # name column
	print "0\t$strand\t"; # score column, strand column 
	print "$otherChrom\t$otherStart\t$otherEnd\t";
	print $otherEnd - $otherStart, "\t"; # otherSize column
	print "$rowNumber\t"; # uid column
	
	print "1000\t"; # posBasesHit
	print 'N/A', "\t", 'N/A', "\t", 'N/A', "\t", 'N/A', "\t"; # testResult, verdict, chits, ccov columns
    print $refcolumnValues->[$colname2Pos{'FILE'}], "\t";
    print $refcolumnValues->[$colname2Pos{'END'}], "\t";
    print $refcolumnValues->[$colname2Pos{'indel_N'}], "\t";
    print $refcolumnValues->[$colname2Pos{'indel_S'}], "\t";
    print $refcolumnValues->[$colname2Pos{'base_S'}], "\t";
    print $refcolumnValues->[$colname2Pos{'base_Match'}], "\t";
    print $refcolumnValues->[$colname2Pos{'base_Mis'}], "\t";
    print $refcolumnValues->[$colname2Pos{'transitions'}], "\t";
    print $refcolumnValues->[$colname2Pos{'transversions'}], "\t";
    print $refcolumnValues->[$colname2Pos{'per_sim'}], "\t";
    print $refcolumnValues->[$colname2Pos{'per_sim_indel'}], "\t";
    print $refcolumnValues->[$colname2Pos{'K_jc'}], "\t";
    print $refcolumnValues->[$colname2Pos{'k_kimura'}], "\n";
			
}

