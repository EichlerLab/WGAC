#!/usr/bin/perl -w

use strict;
use Getopt::Std;


#---------------------------------------------------------------------------------
# the start and end of input file should be sorted with start asc, end desc
#---------------------------------------------------------------------------------

my %opts;
#my $binPrefix   = $ENV{'BINPATH'};
#$binPrefix      .= '/' if ($binPrefix !~ /\/$/ );


if ( !defined $ARGV[0]) {
  print "usage: coordsMerger.pl
        -i the file that has the coordinates to be checked, tab-delimited
        -h if the file for -i has header line

        -s use spaces instead of tab for delimiter
        -u UCSC coordinates convension used (ie, start from 0)

        -n the column for name, ie, merge within the same name group, start column is 0
           if multiple cols are to be used, specify them by comma delimited list. Eg, -n 0,2,3
        -b the column for start coordinate
        -e the column for end coordinate
        -m the keep column other than coordinates columns.
           So merged row will have this value as ':' delimed -m columns from the original rows. Staring from 0.
        -o output file to output merged coordinates
        NOTE: The start and end of input file should be sorted with start asc, end desc
              This merger works for celera coordinates system (the iterator like one)
              so 1080-1089, 1090-1099 will not be merged into 1080-1099\n\n";

  exit;
}


getopts('ushi:n:b:e:o:m:', \%opts);
die "please give a valid column number\n"
  if( !defined $opts{'n'} ||
      !defined $opts{'b'} || $opts{'b'} !~ /^\d+$/ ||
      !defined $opts{'e'} || $opts{'e'} !~ /^\d+$/ ||
      (defined $opts{'m'} && $opts{'m'} !~ /^\d+$/) );


my @kcols = split(/,/, $opts{'n'});
my @orderKey = ();   # to stored keys as their original order
my %all_data = ();
open(IN,  "<$opts{'i'}") || die "can not open input file $opts{'i'}: $!\n";
<IN> if($opts{'h'});
while( <IN> )
  {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;

    next if(length($_) == 0);
    my @cols = ( $opts{'s'} ) ? split(/\s+/) : split(/\t/);


    next if( $cols[$opts{'b'}] !~ /^\d+$/ || $cols[$opts{'e'}] !~ /^\d+$/ );
    my $key = '';
    foreach my $kc (@kcols)
      {
	$key .= $cols[$kc].'@';
      }


    $cols[$opts{'b'}] = $cols[$opts{'b'}] -1 if( !$opts{'u' } );
    #print $cols[$opts{'n'}], ', ', $cols[$opts{'b'}], ', ', $cols[$opts{'e'}], "\n";


    my @addARRAY = ($cols[$opts{'b'}], $cols[$opts{'e'}]);
    push @addARRAY, $cols[$opts{'m'}]  if( defined $opts{'m'} );


    push @{$all_data{ $key } }, [@addARRAY]
      if(defined $all_data{ $key } );

    if(! defined $all_data{ $key } )
      {
	$all_data{ $key } = [ [@addARRAY] ];
	push @orderKey, $key;
      }
  }
close(IN);






open(OUT, ">$opts{'o'}") || die "can not open output file $opts{'o'}: $!\n";
foreach my $key (@orderKey)
  {
    my @sort_data = sort { $a->[0] <=> $b->[0] || $b->[1] <=> $a->[1] } @{$all_data{$key} };
    my ($lastStart, $lastEnd, $lastV, $start, $end, $val);


    foreach my $roA (@sort_data)
      {
	($start, $end) = ($roA->[0], $roA->[1]);
	$val           = $roA->[2] if(defined $opts{'m'} );
	
	# 1st row
	if( !defined $lastStart )
	  {
	    ($lastStart, $lastEnd) = ($roA->[0], $roA->[1]);
	    $lastV                 = $roA->[2] if(defined $opts{'m'} );
	    next;
	  }


	# accual merge within the same name
	if( $start <= $lastEnd )
	  {
	    $lastEnd = $end if( $end > $lastEnd );
	    $lastV   = $lastV.':'.$val if( defined $opts{'m'} && $lastV !~ /\Q$val\E/ );
	    next;
	  }
	
	# $start > $lastEnd or name changed
	if(defined $lastStart)
	  {
	    my @keyVAL = split(/@/, $key);
	    foreach my $kv (@keyVAL)
	      {
		print OUT $kv, "\t";
	      }
	    print OUT $lastStart, "\t", $lastEnd;
	    print OUT "\t", $lastV      if(defined $opts{'m'} );
	    print OUT "\n";
	  }

	($lastStart, $lastEnd) = ($roA->[0], $roA->[1]);
	$lastV                 = $roA->[2] if(defined $opts{'m'} );
      }


    # after last row, which may be merged (not outputting its chunk coords)
    # or merged (its chunk coords not outputted either)
    if(defined $lastStart)
      {
	my @keyVAL = split(/@/, $key);
	foreach my $kv (@keyVAL)
	  {
	    print OUT $kv, "\t";
	  }
	print OUT $lastStart, "\t", $lastEnd;
	print OUT "\t", $lastV      if(defined $opts{'m'} );
	print OUT "\n";
      }

  }
close(OUT);



