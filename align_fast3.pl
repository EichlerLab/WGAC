#!/usr/bin/env perl
####THIS PROGRAM CREATES ALIGNMENTS FROM PARTIAL ALIGNMENTS USING ALIGN PROGRAM  
#ALLOWS FOR LARGE GLOBAL ALIGNMENTS CREATEED QUICKLY####


#LOAD MODULES
use Getopt::Std;
use FindBin;
use strict 'vars';

use vars qw($true $false);
$true=1; $false=0;

use vars qw($opt_i $opt_j $opt_l $opt_o $opt_n $opt_b $opt_w $opt_f $opt_g);
use vars qw($title $position $subseq $header $ppath);
use vars qw($align_path $path1 $path2 $b1 $b2 $e1 $e2);


use vars qw($pusage $program $pversion $pdescription $pgenerate);
### program stats ###
$program = "align_fast";
$pversion='2.010109';
$pdescription = "$program (ver:$pversion)  creates large global alignments quickly from smaller global alignments using the fasta2 program align.";
$pgenerate= 'jeff:dnhc genetics:dnh';
$pusage ="$program -i [file1:b1:e1] -j[file2:b2:e2] [options]";
### program stats end ###

use vars qw($default_gapopening $default_gapextension $default_outputwidth $default_sublength);
$default_gapopening=-16;
$default_gapextension=-4;
$default_outputwidth=60;
$default_sublength=60000;
if ($ARGV[0] eq '') {

print "USAGE
$pusage
DESCRIPTION
$pdescription
ARGUMENTS
-i [path]:[begin]:[end]  path to fasta1:begin(opt):end(opt)
-j [path]:[begin]:[end]  path to fasta2:begin(opt):end(opt)
-l [integer] length subalignments 

   (default $default_sublength bp)
-n [switch] turns off relative positions (not implemented yet)
-f [neg. integer] gap opening penalty
-g [neg. integer] gap extension penalty
-w [integer] bases per line of alignment output 
   (default $default_outputwidth)
-b [integer] fragment bp size if alignments are begin done on fractionated sequence
-o [path] output alignment file (default align.out)
";
exit();
	
}
$ppath="";  #unix bin path
$ppath="perl /binp/" if $^O eq "MSWin32"; #windows bin path
$align_path="$FindBin::Bin/align"; 
getopts('t:i:j:no:f:g:l:b:w:');


$opt_l ||= $default_sublength;
$opt_w ||= $default_outputwidth;
$opt_b ||= '';
if ($opt_b) {
	$opt_b = "-s $opt_b";
}
$opt_i || die "Please enter with -i fasta 1 and :begin:end if necessary.\n";
$opt_j || die "Please enter with -j fasta 2 and :begin:end if necessary.\n";
$opt_f ||= $default_gapopening;
$opt_g ||= $default_gapextension;
print "OUTFILE $opt_o\n";
$opt_o || ($opt_o= 'align.out');
print "########ALIGNSIZE $opt_l\n";
($path1,$b1,$e1)=split ":",$opt_i;
($path2,$b2,$e2)=split ":",$opt_j;
#########CHECK FOR BEGINS AND ENDS REPLACE WITH DEFAULTS IF NOT THERE ##
if ($b1 eq '' || $e1 eq '') {
	
	$e1 =&fasta_length($path1) if $e1 eq '';
	$b1 = 1 if $b1 eq '';
	print "$path1  $b1 $e1\n";
}
if ($b2 eq '' || $e2 eq '') {
	$e2 =&fasta_length($path2) if $e2 eq '';
	$b2 = 1 if $b2 eq '';
	print "$path2 $b2 $e2\n";
}
my ($filename1,$filename2)=($path1,$path2);	
$filename1=~s/^.*\///;
$filename2=~s/^.*\///;
my $orient2='F'; $orient2='R' if $b2>$e2;

	
###########ALIGN THEM###########
my ($alignseq1,$alignseq2)= &align_fast( '-seqpath1'=>$path1, -g=> $opt_g, -f=>$opt_f,'-b1'=>$b1, '-e1'=>$e1,
				'-seqpath2' => $path2, '-b2'=>$b2, '-e2'=> $e2,-frag_size=>$opt_l);
	
#########PRINT THEM ####################
open (OUT,">$opt_o") || die "Can't open outfile ($opt_o)!\n";
&align_save(-filehandle=>\*OUT,-seq1=>$alignseq1,-seq2=>$alignseq2,-name1=>$filename1,
	-name2=>$filename2,-begin1=>$b1,-end1=>$e1,-begin2=>$b2,-end2=>$e2,-orient2=>$orient2, -width=>$opt_w);
close OUT;			

print "$path2\n";
sub align_fast { 		
	my %args= (-seqpath1=>'', -b1=>'',-e1=>'', -seqpath2=>'',
					-b2=>'',-e2=>'', -frag_size=>60000,
					-m1=>700, -m2=>100, -g=>-16, -f =>-4, @_);
	my $r=int(rand(10000000));
	my $first_check=1;
	my $s1_path=$args{'-seqpath1'};
	my $s2_path=$args{'-seqpath2'};
	print "$args{'-seqpath1'}X$args{'-seqpath2'}XS  \n";
	my ($b1,$b2,$e1,$e2) =( $args{'-b1'},$args{'-b2'},$args{'-e1'},$args{'-e2'});
	print "$b1 X $b2 X $e1 X $e2\n";
	my $frag_size=$args{'-frag_size'};
	
	die "No sequence path -seqpath1 ($s1_path) in subroutine align_fast\n" if $s1_path eq '';
	die "No sequence path -seqpath2 ($s2_path) in subroutine align_fast\n" if $s2_path eq '';
	###initial calculations ###
	my $single_piece=$false;
	$single_piece= $true if ($e1-$b1+1 < $frag_size && abs($e2-$b2)+1 <$frag_size );
	my $orient2='F';
	$orient2='R' if $b2>$e2;
	my ($finished1,$finished2);
	my ($align_frag1,$align_frag2,$align_frag1_last,$align_frag2_last);
	my ($mpos0,$mpos1,$mpos2,$mpos2to1)=(0,0,0,0);
	my ($mpos1_last,$mpos2_last,$mpos2to1_last);
	my (@p1,@p2,@p1_last,@p2_last);
	my ($s1_b,$s2_b)=($b1,$b2);
	my ($s1_b_last,$s2_b_last)=(0,0);
	my ($s2_e,$s1_e);
	MAIN:while(1==1) {
		$s1_e=$s1_b+$frag_size-1;
	 	$s1_e=$e1 if $s1_e>$e1;
		if($s1_e>$e1-1000) { $s1_e=$e1; $s2_e=$e2;} #take on a thousand if it can finish
	 	if ($b2<$e2) {
			$s2_e=$s2_b+$frag_size-1;
			$s2_e =$e2 if $s2_e>$e2;
		} else {
			$s2_e=$s2_b-$frag_size+1;
			$s2_e=$e2 if $s2_e<$e2;
	 	}
		#####EXTRACT SUBSEQUENCES AND RUN #######
		my $command="$FindBin::Bin/fasta_subseq33.pl -f $s1_path -o frag1_$r -b $s1_b -e $s1_e -F $opt_b";
		print "SUBSEQ1\n$command\n";
		system "$command";
		$command='';
		$command= "$FindBin::Bin/fasta_subseq33.pl -f $s2_path -o frag2_$r -b $s2_b -e $s2_e -F $opt_b";
		$command = "$FindBin::Bin/fasta_subseq33.pl -f $s2_path -o frag2_$r -b $s2_e -e $s2_b -F -r $opt_b" if $b2>$e2;
		print "SUBSEQ2\n$command\n";
		system $command;
		print "ALIGN\n";
		system "$align_path -f $args{'-f'} -g $args{'-g'}  frag1_$r frag2_$r  > alignout_$r";
		
		($align_frag1_last,$align_frag2_last)=($align_frag1,$align_frag2);
		($align_frag1,$align_frag2)=&align_load(-alignment=> "alignout_$r");
		@p1_last=@p1;
		@p2_last=@p2;
		@p1 = align_position_array(-alignment=>$align_frag1,-start=>$s1_b,-orient=>'F');
		@p2 = align_position_array(-alignment=>$align_frag2,-start=>$s2_b,-orient=>$orient2);
		
		
		####
		my ($max_position,$max_score, $score)=(0,0,0);
		my ($s,$aligned_bases)=(0,0);
		my $gap=$false;
		my $m2_checked=$false;
		for (my $i=0; $i<length($align_frag1); $i++) {
			($s,$gap) = &score(  substr($align_frag1,$i,1),substr($align_frag2,$i,1),$gap );
			#print "$s ";
			 $score+=$s;
			 $aligned_bases++ if $s==1;
			 if ($i==0 &&  $aligned_bases==1 && $mpos0 != 0) {   #postion of mpos2
				print "MPOS2:$p1[$i]==$p1_last[$mpos2] $p2[$i] == $p2_last[$mpos2] MATCHES:$aligned_bases ";
			 	if ($p1[$i]==$p1_last[$mpos2] &&	$p2[$i] == $p2_last[$mpos2]) {
			 		print "AOK\n";
			 	} else {
					die "BAD MPOS2\n";
			 	}
			 }
			 
			 if ($aligned_bases== ($mpos2to1+1) && $mpos0!=0 && $m2_checked==$false) {
			 	$first_check=0; $m2_checked=$true;
			 	print "M2:$args{'-m2'}  ($mpos2to1)\n";
				print substr($align_frag1,0,$i+1),"  ",substr($align_frag1,$i,1),"\n";
				print substr($align_frag2,0,$i+1),"  ",substr($align_frag2,$i,1),"\n";
				print "MPOS1:$p1[$i]==$p1_last[$mpos1] &&	$p2[$i]==$p2_last[$mpos1] MATCHES:$aligned_bases ";
				if($p1[$i]==$p1_last[$mpos1] && $p2[$i]==$p2_last[$mpos1]) {
			 		print "AOK\n";
			 	} else {
	  				print "BAD MPOS1 ALIGNMENT--GOING BACK ONE FRAG\n";
					die "BAD MPOS2\n";			
			 	}
			 }
			if($score>$max_score) {
			  $max_position=$i;
			  $max_score=$score;
			 }
		}
		$mpos0=$max_position;  #this is the first match of the good end#	
		if ($mpos0==0 && !$single_piece) {
			die "BAD ALIGNMENT $mpos0 $max_position $score $max_score\n";
		}
		$mpos1_last=$mpos1; $mpos2_last=$mpos2; $mpos2to1_last=$mpos2to1;
		######END IF AT END OF ALIGNEMN###################
		if ($s1_e==$e1 && $s2_e==$e2) {
			$finished1.=substr($align_frag1_last, 0,$mpos2_last);
			$finished2.=substr($align_frag2_last, 0,$mpos2_last);
			$finished1.=$align_frag1;
			$finished2.=$align_frag2;
			print "LEN:(",length($finished1),") (",length($finished2),")\n";
			unlink ("frag1_$r", "frag2_$r" ,"alignout_$r" );
			return $finished1, $finished2;
		}
		################################################
		##########FIND MPOS1 and MPOS2 from right to left
		
		$mpos1=0; $mpos2=0;
		$gap=$false;
		$aligned_bases=0;
		my ($mpos1_found,$mpos2_found)=($false, $false);
		my $matches_in_a_row=0;
		for (my $i=$mpos0; $i>=0; $i--) {
			($s,$gap) = &score(  substr($align_frag1,$i,1),substr($align_frag2,$i,1),$gap );
			 $score+=$s;
			 $aligned_bases++ if $s==1;
			 if ( ($aligned_bases==$args{'-m1'}) && $s==1 && $mpos1_found==$false) {
			 		 $mpos1=$i;
					 $aligned_bases=0;
					 $mpos1_found=$true;
			 }
			 if (($aligned_bases>= $args{'-m2'}) && $mpos2_found==$false &&$mpos1_found==$true) {
			 	if ($s==1) {$matches_in_a_row++;} else {$matches_in_a_row=0;}
				print "$i: $s $matches_in_a_row\n";
				if ($matches_in_a_row >= 5) {
					$mpos2=$i;
					$mpos2to1=$aligned_bases;
					$mpos2_found=$true;
				}
			 }
						
	  	}
	  if ($mpos1==0 || $mpos2==0 ) {
	  		print "BAD ALIGNMENT--MPOS1 OR MPOS2 equalled zero\n";
			print "INCREASING SEARCH SIZE\n";
	  		@p1=@p1_last; @p2=@p2_last;
			$align_frag1=$align_frag1_last; $align_frag2=$align_frag2_last;
			$mpos1=$mpos1_last; $mpos2=$mpos2_last;
			$mpos2to1=$mpos2to1_last;
			$frag_size*=2;
			$frag_size=50000 if $frag_size>50000;
			next;
		}
		
		###################################################
		##### LAST FRAGMENT #######
		print substr($finished1, length($finished1)-10)," "
				,substr($align_frag1_last,0,10),"\n";
		print substr($finished2, length($finished2)-10), " ",
				substr($align_frag2_last,0,10),"\n";
		$finished1.=substr($align_frag1_last, 0,$mpos2_last);
		$finished2.=substr($align_frag2_last, 0,$mpos2_last);
		#####################
	  	print "ZEROS: $p1[0]  $p2[0]\n";
		print "MPOS 0: $mpos0  1:$mpos1 2:$mpos2\n";
		print "POSs1 0: $p1[$mpos0] 1:$p1[$mpos1] 2:$p1[$mpos2]\n";
		print "POSs2 0: $p2[$mpos0] 1:$p2[$mpos1] 2:$p2[$mpos2]\n";
	  	#################################################
	  	#######SETUP FOR NEXT SEARCH
		$s1_b_last=$s1_b; $s2_b_last=$s2_b; 
		print substr($align_frag1,$mpos2,$mpos1-$mpos2+1),"  ",substr($align_frag1,$mpos1,1),"\n";
		print substr($align_frag2,$mpos2,$mpos1-$mpos2+1),"  ",substr($align_frag2,$mpos1,1),"\n";
	  	$s1_b=$p1[$mpos2];
	  	if ($b2<$e2) {$s2_b=$p2[$mpos2];} else {$s2_b=$p2[$mpos2];}
		print "NEW BEGIN s1: $s1_b   NEWBEGIN S2: $s2_b\n";
		$frag_size=$args{'-frag_size'};
		#print "PAUSE\n";	my $pause=<STDIN>;
		
	  	
	}

}	
sub position_check {
	my ($qb,$qe,$sb,$se,$pc)=@_;
	my @c=@{$pc};
	my $s_rev=$false;
	$s_rev=$true if $sb > $se;

	 return $qb,$qe,$sb,$se;

}
	
sub score {
	my ($b1,$b2,$gap)=@_;
	my $score=0;
	if ($b1 ne '-' && $b2 ne '-') {
		#bases
		$gap=$false;
		if ($b1 eq 'N' || $b2 ne 'N') {
			if ($b1 eq $b2) { $score= 1;} else {$score= -2;}
		} 					
	} else {
		#indel
		$score= -5 if $gap==$false;
		$gap=$true;
	}
	return $score, $gap;
}



sub align_position_array {
	my %args =(-alignment=>'',-start=>1,-orient=>'F',@_);
	my $s1=$args{'-alignment'};
	my $position=$args{'-start'};
	my $orient=$args{'-orient'};
	print "ORIENT $orient\n";
	if ($orient eq 'F'){$position--} else {$position++}
	my @pos_array;
	for (my $i=0; $i<length($s1); $i++) {
		if ( substr($s1,$i,1) ne '-') {
			if ($orient eq 'F'){$position++} else {$position--}
		}
		$pos_array[$i] = $position;
		$pos_array[$i] = '-' if substr($s1,$i,1) eq '-';
	}
	return @pos_array;
}
		

sub align_load { 	
	my %args=(-alignment=>'',@_);
	my ($line, $s1,$s2);
	open (IN,$args{'-alignment'}) || die ("No alignment file ($args{'-alignment'})!");
	$line =<IN> until $line=~/(Global) alignment score:/  || eof(IN);
	if ($1 eq "Global") {
		$line=<IN>; $line=<IN>; $line=<IN>;
		while (!eof(IN) ) {
			$line=~s/\r//;
			chomp $line;
			$line=~ / ([MBDHVRYKSWACTGN-]+)$/;
			$s1.=$1;
			$line=<IN>; $line=<IN>;
			$line=~ / ([MBDHVRYKSWACTGN-]+)$/;
			$s2.=$1;
			$line=<IN>; $line=<IN>; $line=<IN>;	 $line=<IN>;
		}
	} else {
		die "Couldn't load Alignment\n";
	}
	$s1=~s/ +//mg;
	$s2=~s/ +//mg;
	close IN;
	return $s1,$s2;
}
sub align_save {
	my %args=(-filehandle=> \*STDOUT, -seq1=>'seq1', -seq2=>'seq2', -begin1=>1, -end1=>'',-begin2=>1,-end2=>'',
				-name1=>'', -name2=>'', -orient2=>'F', ,-txt=>"-f $opt_f -g $opt_g", -width=>60,@_);
	my $filehandle=$args{'-filehandle'};
	my ($seq1,$b1,$name1,$seq2,$b2,$name2,$orient2,$e1,$e2);
	($seq1,$b1,$name1,$e1)=($args{'-seq1'},$args{'-begin1'},$args{'-name1'}, $args{'-end1'});
	($seq2,$b2,$name2,$orient2,$e2)
			=($args{'-seq2'},$args{'-begin2'},$args{'-name2'},$args{'-orient2'},$args{'-end2'});
	my $width=$args{'-width'};
	my ($pos1,$pos2)=($b1,$b2);
	my ($pos1_del, $pos2_del)=($false,$false);
	#####HEADER ##########
	print $filehandle "FAST ALIGN   $name1 ($b1 to $e1) vs $name2 ($b2 to $e2)\n";
	print $filehandle "Global alignment with ",length($seq1)," spaces. $args{'-txt'}\n";
	print " ";
	####ALIGNMENT ########
	my $m=0;
	for ($m=0; $m+$width<length ($seq1); $m+=$width) {
		
		my $ss1=substr($seq1, $m, $width);
		my $ss2=substr($seq2, $m, $width);
		my $filler='';
		my ($num1,$num2)=('','');
		for (my $i=0; $i<length($ss1); $i++) {
			my ($b1,$b2)=(substr($ss1,$i,1),substr($ss2,$i,1) );
			$pos1++ if $b1 ne '-';
			$pos1_del=$false if $b1 ne '-';
			if ($b2 ne '-') {
				$pos2_del=$false;
				if ($orient2 eq 'F') { $pos2++;} else {$pos2--;}
			}
			if ($pos1 % 10 ==0 && $pos1_del==$false) {
				$pos1_del=$true;
				my $add=' 'x($width+20) . $pos1;
				my $add_len= length ($num1) -$i-11;
				$num1.=substr($add,$add_len);
			}
			if ($pos2 % 10 ==0 && $pos2_del==$false) {
				$pos2_del=$true;
				my $add=' 'x($width+20) . $pos2;
				my $add_len= length ($num2) -$i-11;
				$num2.=substr($add,$add_len);
			}
			if ( $b1 eq  $b2 ) {
				if ($b1 ne 'N' && $b2 ne 'N') {$filler.='|'; } else { $filler.=' ';}
			} else {
				if ( $b1 eq '-' || $b2 eq '-' ) 	{$filler.=' ';} else {	$filler.='*';}
			}
		}
		print $filehandle "\n";
		print $filehandle " $num1\n";
		print $filehandle substr($name1.' 'x10,0,9)," ",$ss1,"\n";
		print $filehandle substr( $m+1 .' 'x10,0 ,10), $filler,"\n";
		print $filehandle substr($name2.' 'x10,0,9)," ", $ss2,"\n";
		print $filehandle " $num2\n";

	}
	if ($m <=length($seq1)) {
		my $ss1=substr($seq1, $m);
		my $ss2=substr($seq2, $m);
		my $filler='';
		my ($num1,$num2)=('','');
		for (my $i=0; $i<length($ss1); $i++) {
			my ($b1,$b2)=(substr($ss1,$i,1),substr($ss2,$i,1) );
			$pos1++ if $b1 ne '-';
			$pos1_del=$false if $b1 ne '-';
			if ($b2 ne '-') {
				$pos2_del=$false;
				if ($orient2 eq 'F') { $pos2++;} else {$pos2--;}
			}
			if ($pos1 % 10 == 0 && $pos1_del == $false) {
				$pos1_del=$true;
				my $add=' 'x ($width+20) . $pos1;
				my $add_len= length ($num1) -$i-11;
				$num1.=substr($add,$add_len);
			}
			if ($pos2 % 10 == 0 && $pos2_del == $false) {
				$pos2_del=$true;
				my $add=' 'x ($width+20) . $pos2;
				my $add_len= length ($num2) -$i-11;
				$num2.=substr($add,$add_len);
			}
			if ( $b1 eq  $b2 ) {
				if ($b1 ne 'N' && $b2 ne 'N') {$filler.='|'; } else { $filler.=' ';}
			} else {
				if ( $b1 eq '-' || $b2 eq '-' ) 	{$filler.=' ';} else {	$filler.='*';}
			}
		}		
		print $filehandle "\n";
		print $filehandle " $num1\n";
		print $filehandle substr($name1.' 'x10,0,9)," ",$ss1,"\n";
		print $filehandle substr( $m+1 .' 'x10,0 ,10), $filler,"\n";
		print $filehandle substr($name2.' 'x10,0,9)," ", $ss2,"\n";
		print $filehandle " $num2\n";
	}
}
sub fasta_length {
	open (FASTA, "$_[0]") || die "File ($_[0]) could not be opened!\n";
	my $header=<FASTA>;
	die ("File does not contain a fasta header!\n") if $header !~/^>/;
	my $totbases=0;
	while ( <FASTA> ) {
		s/\r\n/\n/;
		s/ +//;
		chomp;
		chomp;
		$totbases += length;
	}
	close FASTA;
	return $totbases;
}

