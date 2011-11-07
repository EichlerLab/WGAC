

#!/usr/bin/perl
#
# this is used to sort hte _all table
# after the cat singletone into the table
# sort based on tName[13],and tStart[15]

use strict;

my @a;
my %chrHash;    # this is the hash table based on chr key
my $file = $ARGV[0];
my $firstSortCol = $ARGV[1];
my $secondSortCol = $ARGV[2];
open(F, "$file");
while(<F>){
    if (/query/||/begin/){# header 
		     print $_;
		     next;
    }else{
	# the data
	# collect data based on chr number first
        @a = split(/\s+/, $_);
#print "$a[$firstSortCol]--\n";
        if ($chrHash{$a[$firstSortCol]}){
	    $chrHash{$a[$firstSortCol]} .= $_;
	}else{
		$chrHash{$a[$firstSortCol]} = $_;	
        }
    }
}

sortNprint(\%chrHash);

#------- subs ---------

sub   sortNprint{
# this is a sub to sort the record first and then print them out
#
    my $hashRef = shift;
    my @chr = keys(%$hashRef);
    my $c;
    my ($s,@singleChr, %singleChr, @line, @start, $st);
    
    ## sort the array @chr
    @chr = sort (@chr);
   
    foreach $c(@chr){

#	print "$c-----------\n";
				@singleChr = split(/\n/, $$hashRef{$c});
				foreach $s(@singleChr){
	    		@line = split(/\s+/, $s);

	## you can change second positin you wan tot sort 
	##
					$st = $line[$secondSortCol];
					
  ##---------------------------------------------##
	    		if ($singleChr{$st}){
						$singleChr{$st} = "$singleChr{$st}\n$s";
	    		}else{
							$singleChr{$st} = $s;
	    		}
			}

        # after collecting all lines for singleChr
	# sort them based on the keys
	    @start = keys(%singleChr);
	    @start = sort{$a <=> $b} @start;
	    
	# now to print the records for sa single chromosomeal information 
	# based on the start in after sorted on start

	    foreach $st(@start){
					print "$singleChr{$st}\n";
	    }
	## reset the hashs
			undef %singleChr;
    }
}

