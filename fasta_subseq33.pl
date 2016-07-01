#!/usr/bin/perl -w

#LOAD MODULES
use Getopt::Std;
use strict 'vars';

#based on fasta_subseq2.pl#

#100106 Tin Louie modified sub fasta_getsubseq_frac to look for files with a 3- or 4-digit extension 
#010118 added ability to retrieve from fractionated files#
#010111 added program stats#
#001019 adding batch capablities#

use vars qw($true $false);
use vars qw($opt_f $opt_b $opt_H $opt_h $opt_e $opt_l $opt_s $opt_a $opt_L $opt_D $opt_o $opt_r $opt_c $opt_w $opt_n $opt_C $opt_l $opt_F $opt_N);
use vars qw(@pieces $subseq $header);

use vars qw($program $pversion $pdescription $pgenerate);
$program = "$0";
$program =~ s/^.*\///;
### program stats ###
$pversion='33.010309';
$pdescription = "$program (ver:$pversion) rapidly extras subsequence(s) from whole or fractionated fasta records";
$pgenerate= 'jeff:dnhc genetics:dn';
### program stats end ###



if (! defined $ARGV[0]) {

print "usage: $program -L or -f [path]  [options]
DESCRIPTION\n$pdescription
ARGUMENTS
***FOR BATCH FILES****
-L file1:file2: with list tab-delimited list of accessions begin end 
     in columns 1st,2nd and 3rd columns
  -h no header in list
-H [column] to find name for fasta_file (default none)
-D dir1:dir2 where  fasta files are located
OR
***FOR SINGLE SEQUENCE AND SUBSEQUENCE***
-f [path] fasta sequence file is required
   two of the following -b, -e, -l options are required:
-b [integer] begin position in sequence
-e [integer] end position in sequence
     -b > -e then r	everse complement generated
-l [integer] length of sequence to extract
-s [integer] size of fractionated files
     (default is no fractionation)
-a [integer] append incremeting number to beginning of file name with  this many zero places
-o [path] for output fasta file 
   (otherwise input filename_begin_end is used
   filename is also used for >header if -F invoked
-r [switch] reverse complement (handles atgcyrkmwsbvhdn)
       note!: -b 100 -e 1 -r  reverse the reverse giving you forward
-C [switch] force to uppercase
-c [switch] force to lowercase
-w [switch] number of bases per row (default is 60)
-F [switch] use file name instead of internal header
-N [switch] add numbers in a column (not a true fasta file)
-n [switch] add numbers but same numbering as original file
-v [switch] verbose output for debugging or for a feelin of accomplishment
";

exit;

}
if ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
	system "perldoc $0\n";
	exit;
}

getopts('f:b:e:L:s:l:o:w:rCcFNa:nD:H:h');
$opt_w ||= 60;
$opt_N = 1 if $opt_n;
$opt_r ||= 0;
$opt_L ||= 0;
#print "OPT_S:$opt_s\n";
#######calculate begin, end and length #####
if ($opt_L) {
	die "Please enter outfile with -o ($opt_o)!\n" if !defined $opt_o;
	die "Please enter fasta directories with -D!\n" if !defined $opt_D;
	my @lists;
	(@lists=split":",$opt_L);
	foreach my $l (@lists) {
		open (LIST, $l) || die "Can't open list($l)!\n";
		my $header=<LIST> if !$opt_h;
		while (<LIST>) {
			next if !/\t/;
			s/\r\n/\n/;
			chomp;
			my @c = split "\t";
			print "$c[0] $c[1] $c[2]\n";
			if ($c[0] eq '' || $c[1] !~ /^\d+$/ || $c[2] !~ /^\d+$/ ) {
				 warn "SKIPPING BAD line ($c[0]) column 0)\n";
				 next;
			}
			$c[0].="_000" if $opt_s;
			my $file=&find_file("$opt_D",$c[0]);
			if ($file eq '') {
				 die "BAD FILE NAME  ($c[0]) column 0)\n";
				 next;
			}
			$file=~s/_000$// if $opt_s;
			$c[0] =$file;
			my %h;
			($h{'f'},$h{'b'},$h{'e'})=@c;
			$h{'n'}=$c[$opt_H] if $opt_H;
			$h{'r'}=$opt_r;
			if ($h{'b'} > $h{'e'}) {
				#reverse the reverse#
				($h{'e'}, $h{'b'}) = ($h{'b'} , $h{'e'}) ;
				if ($h{'r'}) {$h{'r'} =0;} else {$h{'r'}=1;};
			}
			print "$h{'f'} $h{'b'}-$h{'e'} R:$h{'r'}\n";
			push @pieces, \%h;
			#my $pause=<STDIN>;
		}
		close LIST;
	}
} else {
	#assume single file inputted#
	if ($opt_b && $opt_e) {
		if ($opt_b > $opt_e) {
			##reverse it
			$opt_r||=0;
			$opt_r= !$opt_r;
			($opt_b,$opt_e)=($opt_e,$opt_b);
		}
		$opt_l = $opt_e - $opt_b + 1;
	} elsif ($opt_b && $opt_l) {
		$opt_e = $opt_b + $opt_l - 1;
	} elsif ($opt_e && $opt_l) {
		$opt_b = $opt_e - $opt_l + 1; 
	} else {
		die "Must specify at least 2 of begin(-b), end(-e) or length(-l) to extract\n";
	}
	$pieces[0]{'f'}=$opt_f;
	$pieces[0]{'b'}=$opt_b;
	$pieces[0]{'e'}=$opt_e;
	$pieces[0]{'r'}=$opt_r;
}

if ($opt_L && !$opt_o) {
	print "Put files in current directory\n";
	$opt_o= '.';
}
if ($opt_o && $opt_L) {
	if ($opt_L && -d $opt_o) {
	 		###directory to put multiple files###
	} else {
		###file to store all output###
		open(OUT, ">$opt_o") || die "Can't open out file!\n";; 
	}
} else {
	if ($opt_o ) {  #file#
		
	} else {
		$opt_o=$opt_f;
		$opt_o =~ s/^.*\///;
		$opt_o = "./$opt_o.$opt_b.$opt_e";
	}
	open(OUT, ">$opt_o") || die "Can't write to ($opt_o)!\n";;
}
print "OUT:$opt_o\n";
############################
############################
#####MAIN LOOP       #######
my $count;
foreach my $p (@pieces) {
	my ($f,$b,$e,$r)=($$p{'f'},$$p{'b'},$$p{'e'},$$p{'r'});
	print "$f  $b $e $r\n";
	$opt_l=$e-$b+1;
	if ($opt_s) {
		&fasta_getsubseq_frac($f,$b,$e);
	} else {
		&fasta_getsubseq_whole($f,$b,$e);
	}
	my $orig_header=$header;
	if ($r) {
		$subseq = reverse $subseq;
		$subseq =~ tr/atgcyrkmwsbvhdATGCYRKMWSBVHD/tacgrymkswvbdhTACGRYMKSWVBDH/; 
		($b,$e)=($e,$b);
	}
	####change case #########
	$subseq=uc $subseq if $opt_C;
	$subseq=lc $subseq if $opt_c;
	
	my $header="$$p{'f'}.$$p{'b'}.$$p{'e'}";
	$header.="_R" if $r;
	$header =~ s/^.*\///;
	
	if ($opt_F) {
		$orig_header =~ s/>//;
		if ($opt_f) {
			my $head= $opt_o;
			$head=~ s/.*\///;
			$header =$head;
		} else {
			$orig_header =~ s/>//;
			$header.=" $orig_header";
			chomp $header;
		}
	} else {
		$header=$1 if $orig_header =~ />(\S+)/;
		#print "\n$header\n";
		$header=~s/_\d\d\d$/ / if $opt_s;
	}
	if ($opt_a) {
		$count++;
		
		my $append=substr(("0" x $opt_a) ."$count",-$opt_a);
		print "TRYING TO APPEND NUMBER $append\n";
		$header ="$append$header";
	}
	if ($opt_H) {
		$header=$$p{'n'};
	}
	if ($r) {$r='R';} else {$r='F';}
	if ($opt_L && -d "$opt_o" ) {
		close OUT;
		if ($opt_H) {
			open (OUT, ">$opt_o/$$p{'n'}") || die "CAn't create ($opt_o/$$p{'n'})!\n";
		} else {
			$f=~s/^.*\///;
			open (OUT, ">$opt_o/$f"."_$b"."_$e") || die "Can't create ($opt_o/$f"."_$b"."_$e)!\n";
		}
	}
	print ">$header from $b to $e  ($opt_l)($r)\n";
	print OUT ">$header from $b to $e  ($opt_l)($r)\n";	
	if ($opt_n) {
		print OUT &breaksequence($subseq, $opt_w, $opt_N, $opt_b-1);
	} else { 
		print OUT &breaksequence($subseq, $opt_w, $opt_N, 0)
	}

}
close OUT;
################################
################################
####### SUBROUTINES ############
################################

sub fasta_getsubseq_whole {
	my $file=shift;
	my $begin=shift;
	my $end=shift;
	my $position=0;
	$subseq='';
	open (FASTA,"$file") || die "Can not open fasta file($file)!\n";
	$header = <FASTA>;
	$header =~ s/\r\n/\n/;
	chomp $header;
	while ( <FASTA> ) {
		s/\r\n/\n/;
		chomp; 
		s/\s+//;
		$position += length;
		if ($position >= $begin) {
			my $start= $begin - ($position- length ) -1;
			$start=0 if $start<0;
			#print "START:$start\n";
			$subseq.= substr($_,$start);

		}
		last if $position >= $end;
	}
	close FASTA;
	die "Trim ($begin to $end []) too long! Length is $position\n" if $end > $position;
	$subseq = substr($subseq, 0, $end -$begin+1 );
}

#####create header ######
sub find_file {
	my ($paths, $names) = @_;
	my @paths=split ":",$paths;
	my @names=split ":",$names;
	for my $path (@paths) {
	 	for my $name (@names) {
	 		my $p = $path;
	 		$p .= "/$name";
	 		#print "TESTING $p\n";
	 		return $p if (-e $p);
	 	}
	}
	return "";
}




sub fasta_getsubseq_frac {
	my $opt_f=shift;
	my $opt_b=shift;
	my $opt_e=shift;
	my $bfasta=int(($opt_b-1)/$opt_s);
	my $bpos=$opt_b-$bfasta*$opt_s;
	my $efasta=int(($opt_e-1)/$opt_s);
	my $epos=$opt_e-$efasta*$opt_s;
	print "$bfasta:$bpos:$efasta:$epos\n";
	$subseq='';
	$header='';
	for (my $i=$bfasta; $i<=$efasta; $i++) {
		my $path= $opt_f."_" . substr("000$i",-3);

		# if file with 3-digit extension does not exist, then assume file has 4-digit extension
		$path = $opt_f . "_" . substr("000$i",-4) unless (-r $path);

		my $begin=1;
		my $end=$opt_s;
		$begin=$bpos if $i==$bfasta;
		$end=$epos if $i==$efasta;
		if ($begin != 1 || $end != $opt_s) {
			print "EXTRACTING SUBSEQ:$path:$begin:$end(I:$i)\n";
			$subseq.= &getsubseq_frac($path,$begin,$end);
		} else {
			print "EXTRACTING WHOLESEQ:$path:$begin:$end(I:$i)\n";
			$subseq.= &getwholeseq_frac ($path);
		}
		print "SUBSEQ LEN:",length($subseq);
	}
}


sub getsubseq_frac {
	my $fasta=shift;
	my $begin=shift;
	my $end=shift;
	my $subseq='';
	my $position=0;
	#this sub returns nothing because it is built for speed#
	open (FASTA,$fasta) || die "Can not open fasta file $fasta\n";
	$header=<FASTA>;
	die "Improper fasta $fasta\n" if $header !~/^>/;
	while ( <FASTA> ) {
		s/\r\n/\n/;
		chomp; chomp;
		s/\s+//;
		$position += length;
		if ($position >= $begin) {
			my $start= $begin - ($position- length ) -1;
			$start=0 if $start<0;
			#print "START:$start\n";
			$subseq.= substr($_,$start);

		}
		last if $position >= $end;
	}
	die "Extracting $fasta ($begin to $end)failed! Request to long!(end length:$position)!\n" if $end > $position;
	$subseq = substr($subseq, 0, $end -$begin+1 );
	close FASTA;
	return $subseq
}

sub getwholeseq_frac {
	my $position=0;
	my $fasta=shift;
	open (FASTA,$fasta) || die "Can not open fasta file $fasta\n";
	$header=<FASTA>;
	die "Improper fasta $fasta\n" if $header !~ /^>/;
	while ( <FASTA> ) {
		s/\r\n/\n/;
		chomp;
		s/\s+//;
		$position += length;
		$subseq.=$_;
	}
	die "End $position not equal to $opt_s\n!" if $opt_s !=$position
}

sub breaksequence {
	my ($whole_seq, $width, $numb_on, $plus_num) =@_;
	my $m=0;
	my $broken_seq="";
	for ($m=0; $m+$width<length ($whole_seq); $m+=$width) {
		$broken_seq .= substr($whole_seq, $m, $width);
		$broken_seq .= "  ".($m+$width+$plus_num) if ($numb_on);
		$broken_seq .= "\n";
	}
	if ($m <=length($whole_seq)-1) {
		$broken_seq .=  substr($whole_seq, $m) ;
		$broken_seq .= "  ".length($whole_seq+$plus_num) if $numb_on;
		$broken_seq .= "\n";
	}
	return $broken_seq;
}

