#!/bin/env perl


use strict;
use Getopt::Std;
use DBI;



my $dbug = 0;
my %opts;


if ( !defined $ARGV[0]) {
  print "usage: $0:
        outputs the overlap/exclude of two tables using merge sort algorithm

        -b database
        -l login
        -p pswd
        (only needed when data src are from database tables)

        -i the table for 1st portion with chrom, chromStart, chromEnd
        -j the table for 2nd portion with chrom, chromStart, chromEnd
        -f option for -i comes from file with chrom, chromStart, chromEnd as first 3 columns
        -t option for -j comes from file with chrom, chromStart, chromEnd as first 3 columns
        -c condition for table i. Don't set it if no.
           Given as tableName.column=value. Concate with ' AND ' if multiple
        -d condition for table j. Don't set it if no
        -L output the exclude from left (-i) table (part of -i table not overlapped with -j table)
        -R output the exclude from right(-j) table (part of -j table not overlapped with -i table)
           If -L, -R are set simultaneously, only -L is used
        -o output file


        Eg:
        twoPartOvp_mgsrt.pl         -i 1.tab -f             -j  2.tab -t          -o /tmp/both.tab
        twoPartOvp_mgsrt.pl         -i 1.tab -f             -j  2.tab -t     -L   -o /tmp/onlyIn1.tab
        twoPartOvp_mgsrt.pl         -i 1.tab -f             -j  2.tab -t     -R   -o /tmp/onlyIn2.tab

        NOTE: There can be no overlap within tables itself, ie, the tables have to have N(on)R(edundant) bases for it to work\n\n";

  exit;
}


getopts('LRftb:l:p:i:j:o:c:d:', \%opts);
die "please set either L or R or neither of them. Not both\n" if( defined $opts{'L'} &&  defined $opts{'R'});
my $left  = $opts{'L'};
my $right = $opts{'R'};



my %fstlist = ();
my %sndlist = ();
my $iCond   = $opts{'c'};
my $jCond   = $opts{'d'};
my $dbh;


if( !$opts{'f'} || !$opts{'t'} )
  {
    $dbh = DBI->connect("DBI:mysql:$opts{'b'}", $opts{'l'}, $opts{'p'}, {AutoCommit => 1} ) || die "Couldn't connect to database: ".DBI->errstr;
  }


if( !$opts{'f'} )
  {
    my $statement = "SELECT chrom, chromStart, chromEnd FROM $opts{'i'} ";
    $statement   .= " WHERE $iCond " if( defined $iCond);
    $statement   .= ' ORDER BY chrom, chromStart;';
    my $listRoA = $dbh->selectall_arrayref($statement, { Columns=>{} });
    foreach my $listRoH (@$listRoA)
      {
	push @{ $fstlist{$listRoH->{'chrom'}} }, [$listRoH->{'chromStart'}, $listRoH->{'chromEnd'}]
	  if(defined $fstlist{ $listRoH->{'chrom'} } );

	$fstlist{$listRoH->{'chrom'}} = [ [$listRoH->{'chromStart'}, $listRoH->{'chromEnd'}] ]
	  if(!defined $fstlist{ $listRoH->{'chrom'} } );
      }
  }
else
  {
    open(IN, $opts{'i'}) || die "cannot open $opts{'i'}: $!\n";
    while( <IN> )
      {
	chomp;
	my @data = split(/\t/);
	next if( $data[1] !~ /^\d+$/ || $data[2] !~ /^\d+$/ );

	push @{ $fstlist{$data[0]} }, [$data[1], $data[2]]
	  if(defined $fstlist{$data[0]} );

	$fstlist{$data[0]} = [ [$data[1], $data[2]] ]
	  if(!defined $fstlist{$data[0]} );
      }
    close(IN);
  }




if( !$opts{'t'} )
  {
    my $statement = "SELECT chrom, chromStart, chromEnd FROM $opts{'j'} ";
    $statement .= " WHERE $jCond " if( defined $jCond);
    $statement .= ' ORDER BY chrom, chromStart;';
    my $listRoA = $dbh->selectall_arrayref($statement, { Columns=>{} });
    foreach my $listRoH (@$listRoA)
      {
	push @{ $sndlist{$listRoH->{'chrom'}} }, [$listRoH->{'chromStart'}, $listRoH->{'chromEnd'}]
	  if(defined $sndlist{ $listRoH->{'chrom'} } );

	$sndlist{$listRoH->{'chrom'}} = [ [$listRoH->{'chromStart'}, $listRoH->{'chromEnd'}] ]
	  if(!defined $sndlist{ $listRoH->{'chrom'} } );
      }
  }
else
  {
    open(IN, $opts{'j'}) || die "cannot open $opts{'j'}: $!\n";
    while( <IN> )
      {
	chomp;
	my @data = split(/\t/);
	next if( $data[1] !~ /^\d+$/ || $data[2] !~ /^\d+$/ );

	push @{ $sndlist{$data[0]} }, [$data[1], $data[2]]
	  if(defined $sndlist{$data[0]} );

	$sndlist{$data[0]} = [ [$data[1], $data[2]] ]
	  if(!defined $sndlist{$data[0]} );
      }
    close(IN);
  }
$dbh->disconnect() if( !$opts{'f'} || !$opts{'t'} );




my $union = join(' ', sort keys %fstlist);
foreach my $key (sort keys %sndlist)
  {
    next if( $union =~ /\Q $key \E/ );
    next if( $union =~ /^\Q$key\E/ );
    next if( $union =~ /\Q$key\E$/ );
    $union .= ' '.$key;
  }
print $union,"\n" if( $dbug );
my @allChrom = split(/ /, $union);





open(OUT, ">$opts{'o'}") || die "cannot open $opts{'o'} to write: $!\n";
foreach my $chrom (@allChrom)
  {

    @{$fstlist{$chrom}} = sort {$a->[0] <=> $b->[0] } @{$fstlist{$chrom}}
		if defined $fstlist{$chrom};
    @{$sndlist{$chrom}} = sort {$a->[0] <=> $b->[0] } @{$sndlist{$chrom}}
		if defined $sndlist{$chrom};


    # what is in i but not in j
    if($left)
      {
	next if( !defined $fstlist{$chrom} );
	if( !defined $sndlist{$chrom} )
	  {
	    outputList( $fstlist{$chrom}, $chrom, *OUT);
	    next;
	  }
      }


    # what is in j but not in i
    if($right)
      {
	next if( !defined $sndlist{$chrom} );
	if( !defined $fstlist{$chrom} )
	  {
	    outputList( $sndlist{$chrom}, $chrom, *OUT);
	    next;
	  }
      }


    # what is in i and in j
    if( !$left && !$right )
      {
	next if(  !defined $fstlist{$chrom} || !defined $sndlist{$chrom} );
      }



    # chrom is defined for both i and j
    my $i = 0;
    my $j = 0;
    my ($leftStart, $rightStart) = (-1, -1);
    while( $i < scalar( @{ $fstlist{$chrom} } ) && $j < scalar( @{ $sndlist{$chrom} } ))
      {
	my $fstRoA = $fstlist{$chrom}->[$i];
	my $sndRoA = $sndlist{$chrom}->[$j];



	if( $fstRoA->[1] <= $sndRoA->[0] )
	  {
	    print OUT $chrom, "\t", $fstRoA->[0], "\t", $fstRoA->[1], "\n" if($left && $leftStart == -1 && $fstRoA->[1] > $fstRoA->[0] );
	    print OUT $chrom, "\t", $leftStart,   "\t", $fstRoA->[1], "\n" if($left && $leftStart != -1 && $fstRoA->[1] > $leftStart);
	    $leftStart = -1;
	    $i++;
	    next;
	  }
	if( $sndRoA->[1] <= $fstRoA->[0] )
	  {
	    print OUT $chrom, "\t", $sndRoA->[0], "\t", $sndRoA->[1], "\n" if($right && $rightStart == -1 && $sndRoA->[1] > $sndRoA->[0]);
	    print OUT $chrom, "\t", $rightStart,  "\t", $sndRoA->[1], "\n" if($right && $rightStart != -1 && $sndRoA->[1] > $rightStart);
	    $rightStart = -1;
	    $j++;
	    next;
	  }

	# pick the greatest start and least end
	my $start = ($fstRoA->[0] > $sndRoA->[0]) ? $fstRoA->[0] : $sndRoA->[0];
	my $end   = ($fstRoA->[1] > $sndRoA->[1]) ? $sndRoA->[1] : $fstRoA->[1];


	print OUT $chrom, "\t", $start, "\t", $end, "\n" if( !$left && !$right && $end > $start );
	if( $left && $fstRoA->[0] < $start )
	  {
	    print OUT $chrom, "\t", $fstRoA->[0], "\t", $start, "\n" if( $leftStart == -1 && $start > $fstRoA->[0] );
	    print OUT $chrom, "\t", $leftStart,   "\t", $start, "\n" if( $leftStart != -1 && $start > $leftStart );
	  }
	if( $right && $sndRoA->[0] < $start )
	  {
	    print OUT $chrom, "\t", $sndRoA->[0], "\t", $start, "\n" if( $rightStart == -1 && $start > $sndRoA->[0] );
	    print OUT $chrom, "\t", $rightStart,  "\t", $start, "\n" if( $rightStart != -1 && $start > $rightStart  );
	  }


	# move the one whose end is smaller
	$i++ if( $fstRoA->[1] == $end );
	$j++ if( $sndRoA->[1] == $end );
	$leftStart  = ($fstRoA->[1] > $end) ? $end : -1;
	$rightStart = ($sndRoA->[1] > $end) ? $end : -1;
      }

    # append those that are left from -i table
    while($left && $i < scalar( @{ $fstlist{$chrom} } ) )
      {
	if( $leftStart != -1 )
	  {
	    print OUT $chrom, "\t", $leftStart,   "\t", $fstlist{$chrom}->[$i]->[1], "\n" if( $fstlist{$chrom}->[$i]->[1] > $leftStart );
	    $leftStart = -1;
	  }
	else
	  {
	    print OUT $chrom, "\t", $fstlist{$chrom}->[$i]->[0], "\t", $fstlist{$chrom}->[$i]->[1], "\n";
	  }
	$i++;
      }

    # append those that are left from -j table
    while($right && $j < scalar( @{ $sndlist{$chrom} } ) )
      {
	if( $rightStart != -1 )
	  {
	    print OUT $chrom, "\t", $rightStart,   "\t", $sndlist{$chrom}->[$j]->[1], "\n" if( $sndlist{$chrom}->[$j]->[1] > $rightStart );
	    $rightStart = -1;
	  }
	else
	  {
	    print OUT $chrom, "\t", $sndlist{$chrom}->[$j]->[0], "\t", $sndlist{$chrom}->[$j]->[1], "\n";
	  }
	$j++;
      }
  }
close(OUT);



sub outputList{
  my ($ref_array, $chrom, $OUT) = @_;

  foreach my $roA (@$ref_array)
    {
      print $OUT $chrom, "\t", $roA->[0], "\t", $roA->[1], "\n";
    }
}







