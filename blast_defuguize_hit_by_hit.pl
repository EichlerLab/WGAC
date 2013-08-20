#!/usr/bin/perl
###00-06-01  RENAMED FROM blast_parse_hit_by_hit2.pl to blast_defuguize_hit_by_hit.pl
###version two adds the ability to have records in multiple locations/directories


#!/usr/local/bin/perl
#LOAD MODULES

use Getopt::Std;
use strict 'vars';

use vars qw($true $false);
$true=1; $false=0;

use vars qw($opt_o $opt_t $opt_s $opt_d);
use vars qw($title $position $subseq $header $count);

if ($ARGV[0] eq '') {

	print "blast_defuguize_hit_by_hit.pl ******************************************\n";
	print "  This program defuguizes a pairwise blast table.\n";
	print "\t-t [path] a blast parse table file is required\n";
	print "\t-d [path] for the directory(s) with files or encode and mask_out dir \n";
	print "       separate paths with colon\n";
	print "\t-o [path] output table file (default --in.defugu\n";
	print "\t-g [switch] requires only -f and adds .out .fugu\n";
	print "\t-s [switch] do not do subject just query defuguize\n";
	die "***********************************************************\n";
}

getopts('t:o:d:s');
$opt_t || die "Please enter with -t the path for the pairwise blast table\n";
$opt_d || die "Please enter with -d the paths to directory containing  files or directories fugu and mask_out\n";
$opt_o || ($opt_o = $opt_t ."\.defugu");
open (TABLE, $opt_t) || die "Can't open pairwise blast table: $opt_t\n";
open (NEWTABLE, ">$opt_o") || die "Can't output table file: $opt_o\n";
my $table_header=<TABLE>;
print NEWTABLE "$table_header";
$count=0;
while (<TABLE>) {
	print $count++, "\n";
	chomp;
	
	my @c=split "\t";
	#####GET ORIGINAL LENGTH#####
	foreach my $i ( 0,4) {
		next if $i==4 && $opt_s==$true;
		my $path = &find_file($opt_d, 'fugu',$c[$i], '.fugu');
		my ($tmp,$h)=&load_fasta("$path");
		($c[$i+3])=$h=~/UELEN:(\d+)$/;	
		my $switch=$false;
		my ($b,$e)=@c[$i+1,$i+2];
		if ($b > $e) {
			$switch=$true;
			($b,$e)=($e,$b);
		}
		$path=&find_file($opt_d, 'mask_out',$c[$i], '.out');
		($b,$e)=&rm_out_decode(begin => $b,end=>$e,rm_out=>"$path");
		($b,$e)=($e,$b) if $switch;
		@c[$i+1,$i+2]=($b,$e);
	}
	#print "AFTER :\t",join ("\t",@c),"\n";
	print NEWTABLE join("\t",@c),"\n";
}




####decode test ####
#my ($b,$e) = &rm_out_decode(5,30437, $opt_r);
#print "$b to $e \n";


		
sub rm_out_decode {
	my %args=(begin=>'',end=>'',rm_out=>'', @_);
	#cprint "$args{'begin'}  $args{'end'}\n";
	my @repeats = &load_rm_out($args{'rm_out'});
	###need more code to do this
	foreach my $r (@repeats) {
		my $l= $$r{'e'}-$$r{'b'} +1;
		$args{'begin'} += $l if $args{'begin'} >= $$r{'b'} ;
		 if ($args{'end'} >= $$r{'b'}) {
				$args{'end'} += $l;
		} else {
			last;
		}
	}
				
	return ($args{'begin'}, $args{'end'});
}

sub load_fasta {
	open (FASTAIN, "$_[0]") || die "Can't open $_[0]\n";
	my $fasta = '';
	my $header = <FASTAIN>;
	chomp $header;
	while (<FASTAIN>) {
		s/\r\n/\n/;
		chomp;
		$fasta .= $_;
	}
	return ($fasta, $header);
	
}
	
sub load_rm_out {
		open (RMOUT, "$_[0]")  || die "Can't open $_[0]\n";
		my @repeats=();
		my $line =<RMOUT>; $line=<RMOUT>; $line=<RMOUT>;
		my $i;
		while (<RMOUT>) {
			###need to read the file backward
			chomp;
			s/^ +//;
			my @col=split " +";
			if ($i>0 && ($col[5] <= $repeats[$i-1]{'e'} ) ) {
				#print "OVERLAP\n";
				$i--;
			} else {
				$repeats[$i]{'b'} = $col[5];
			}
		
			$repeats[$i]{'e'} = $col[6];
			$repeats[$i]{'rep'} .= "_$col[9]";
			$repeats[$i]{'fam'} .= "_$col[10]";
			#print $repeats[$i]{'b'}, "\n";
			$i++;
		}
		close RMOUT;
	   return @repeats;
}

sub find_file {
	my ($paths, $sub_paths, $names, $extensions) = @_;
	my @paths=split ":",$paths;
	my @sub_paths=split ":",$sub_paths;
	my @names=split ":",$names;
	my @extensions=split ":",$extensions;
	for my $path (@paths) {
	 for my $sub_path (@sub_paths) {
	 	for my $name (@names) {
	 		for my $ext (@extensions) {	
	 			my $p = $path;
	 			$p .= "/$sub_path" if $sub_path;
	 			$p .= "/$name$ext";
	 			print "TESTING $p\n";
	 			return $p if (-e $p);
	 		}
	 	}
	 }
	}
	return "";
}


sub save_fasta {
	my ($filename, $whole_seq, $header) = @_;
	open (FASTAOUT, ">$_[0]") || die "Can't create $_[0]\n";
	print FASTAOUT "$header\n";

	my $width=60;
	my $m=0;
	for ($m=0; $m+$width<length ($whole_seq); $m+=$width) {
		print  FASTAOUT substr($whole_seq, $m, $width),"\n";
	}
	if ($m <length($whole_seq)-1) {
		print FASTAOUT  substr($whole_seq, $m),"\n";
	}
	close FASTAOUT;
}





