#!/usr/bin/perl -w
use strict 'vars';
use Getopt::Std;
use POSIX ":sys_wait_h"; 

$| = 1;
# usage:
# perl fugumation.pl -i fugu -o selfblast
# 	fugu is name of directory containing fugu files
# 	selfblast is name of directory to hold output files

use vars qw($opt_i $opt_o);
getopts('i:o:');

#die "cannot create subdir lav_int \n" unless (mkdir "$opt_o/lav_int");
#die "cannot create subdir lav_int2 \n" unless (mkdir "$opt_o/lav_int2");
#die "cannot create subdir data \n" unless (mkdir "$opt_o/data");

our %aPids;


runWebb(-75, 14, 1400);
runWebb(-75, 14, 1200);
runWebb(-75, 14, 1000);
runWebb(-75, 14, 900);
runWebb(-75, 14, 800);
runWebb(-75, 14, 600);
runWebb(-80, 14, 400);
runWebb(-80, 14, 200);
runWebb(-90, 18, 200);

######BREAK OVERLAPPING LAVS #######
# system "perl blast_lav_break_self_overlap2.pl --in $opt_o/lav_int --out $opt_o/lav_int2";

#####PARSE THE WONDERFUL WORLD OF LAV ########
# system "perl blast_lav_hit_by_hit.pl --in $opt_o/lav_int2 --out $opt_o/data/lav_int2.parse -options 'MIN_BPALIGN=>200, MIN_FRACBPMATCH=>0.88, MAX_%GAP => 40, SKIP_SELF => 0, SKIP_OVERLAP=>1'";

exit;


sub startJob {
    my ( $f, $outputFile, $IVparameter, $Wparameter, $Yparameter ) = @_;


	my $szCommand = "time ./webb_self $opt_i/$f B=2 M=30 I=$IVparameter V=$IVparameter O=180 E=1 W=$Wparameter Y=$Yparameter > $outputFile";

    my $nPid = fork();
    if ( $nPid != 0 ) {
        $aPids{ $nPid } = $szCommand
    }
    else {
        # I am the child process

        print STDERR "DOING $f $szCommand\n";
        system( $szCommand );
        print STDERR "completed $szCommand\n";
        
        # at this point, it has completed or crashed

	    # if output file exists but its first line is "bad", then delete the file
	    # so that it can be redone in the next round

        # debugging
        if ( ! -e $outputFile ) {
            print STDERR "$outputFile doesn't exist\n";
        }
        else {
            print STDERR "$outputFile exists\n";
        }            
        # end debugging

	    if (open (POSTCHECK, $outputFile)) {
	    	my $line= <POSTCHECK>;
	    	close POSTCHECK;
            if ( !defined( $line ) ) {
                print STDERR "no lines in $outputFile\n";
                unlink $outputFile
            }
            else {
                unlink $outputFile unless ( $line =~ /\#/ );
            }
	    }

	    # if webb_self had a core dump, remove the file (it can grow big)
	    unlink 'core' if (-e 'core');



        exit( 0 );
    }
}


sub cpuAvailableOrWait {

    # $# returns the maximum index so 0 means there is 1 element
    my $nLength = keys %aPids;
    if ( $nLength < 2 ) {
        return;
    }

    my $nWaitTimes = 0;

    my $bFoundDeadPid = 0;
    while( !$bFoundDeadPid ) {
        my $nPid;
        foreach $nPid (keys %aPids ) {
            my $nDeadPid = waitpid( $nPid, &WNOHANG );
            if ( $nDeadPid != 0 ) {
                # found dead child process
                # remove it from %aPids


                delete $aPids{ $nDeadPid };
                $bFoundDeadPid = 1;
                last;
            }
        }

        if ( !$bFoundDeadPid ) {
            $nWaitTimes += 1;
            if ( ( $nWaitTimes == 1) || ($nWaitTimes == 2) || ($nWaitTimes == 4) || ($nWaitTimes == 8) || ($nWaitTimes == 16) || ($nWaitTimes == 32) || ($nWaitTimes == 64) || ($nWaitTimes == 128) || ($nWaitTimes == 256) || ($nWaitTimes == 512) || ($nWaitTimes == 1024) || ($nWaitTimes == 2048) || ( $nWaitTimes % 2048 == 0 ) ) {
                print STDERR "waiting... $nWaitTimes\n";
            }
            sleep 1;
        }
    }
    

}


sub waitForLastProcessToFinish {

    my $nLength = keys %aPids;
    print STDERR "waitForLastProcessToFinish length = " . $nLength . "\n";

    my $nPid;

    # debugging
    foreach $nPid (keys %aPids ) {
        print STDERR "still has " . $nPid . "\n";
    }
    # end debugging

    $nLength = keys %aPids;
    if ( $nLength == 0 ) {
        return;
    }

    foreach $nPid (keys %aPids ) {
        $nLength = keys %aPids;
        print STDERR "waiting for pid " . $nPid . " length = " . $nLength . "\n";
        my $nDeadPid = waitpid( $nPid, 0 );
        delete $aPids{ $nPid };
        $nLength = keys %aPids;
        print STDERR "returned from waitpid for " . $nPid . " now length = " . $nLength . "\n"
    }

}







# webb_batch.pl was boiled down to the following subroutine:
# invokes the webb_self program on each file in fugu directory, hopefully creates corresponding output file
# ... won't create output file if it already exists (from a previous iteration)
# ... removes 'core' file which is the result of segmentation fault
sub runWebb {
	my $IVparameter = shift; # penalty
	my $Wparameter = shift; # word size
	my $Yparameter = shift; 

    # notice that we wait for all pids to finish before the next runWebb
    # when we start runWebb, we start with no pids
    our %aPids = ();


	opendir (DIR, $opt_i) || die "cannot open directory to fugu files\n";

	while ( defined (my $f = readdir(DIR)) ) {
	    next if $f =~ /^\./;
		my $outputFile = "$opt_o/$f.intf";
		
		next if (-e $outputFile); # OUTPUT FILE ALREADY EXISTS

        startJob( $f, $outputFile, $IVparameter, $Wparameter, $Yparameter );

        cpuAvailableOrWait();
	}

	closedir DIR;

    waitForLastProcessToFinish();
	
}
